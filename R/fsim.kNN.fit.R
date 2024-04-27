fsim.kNN.fit<- function(x,y,
seed.coeff=c(-1,0,1), order.Bspline=3, nknot.theta=3, knearest=NULL,
min.knn=2, max.knn=NULL, step=NULL,  
kind.of.kernel="quad",range.grid=NULL,nknot=NULL,n.core=NULL)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
if (is.null(max.knn)) max.knn <- (n %/% 2)
p <- ncol(x)
estimated.Y <- list()
length.curve.y<-ncol(y)
yhat.cv <- matrix(0,n,length.curve.y)
if (is.null(range.grid)) range.grid <- c(1,p)
t0 <- mean(range.grid)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
if (is.null(max.knn)) max.knn <- n%/%2
if (is.null(knearest)){
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
k.seq <- knearest
kmax <- max(k.seq)       
num.band <- length(k.seq)
dim.base.theta <- order.Bspline + nknot.theta
THETA.seq <- permutations(length(seed.coeff), dim.base.theta, seed.coeff, repeats.allowed=TRUE)
THETA.seq <- THETA.seq[(rowSums(abs(THETA.seq)) != 0)  & (!is.na(rowSums(abs(THETA.seq)))), ]
a <- range.grid[1]
b <- range.grid[2]
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))     
point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142)
weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
x.gauss<- 0.5 * ((b + a) + (b - a) * point.gauss)
lx.gauss<- length(x.gauss)
Bspline.g.theta<-splineDesign(delta.theta, x.gauss, order.Bspline)
H.theta <-t(Bspline.g.theta)%*%(Bspline.g.theta*(weight.gauss*0.5*(b-a)))
x0 <- seq(a, b, length = p)
Knot <-seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]  
delta <-sort(c(rep(c(a, b),order.Bspline), Knot))
Bspline <-splineDesign(delta,x0,order.Bspline)
Cmat <-crossprod(Bspline)
Dmat1 <-crossprod(Bspline, t(x))
coef.mat1 <-symsolve(Cmat, Dmat1)
Bspline.theta<-splineDesign(delta.theta,x0,order.Bspline)
Bspline.g<-splineDesign(delta, x.gauss, order.Bspline)
H.x <-t(Bspline.g)%*%(Bspline.g*(weight.gauss*0.5*(b-a)))
dim<-ncol(THETA.seq)
per<-nrow(THETA.seq)
theta.normal<-matrix(0,nrow=per,ncol=dim)
for(i in 1:per){
	coefi=THETA.seq[i,]
	prod.coef <-outer(coefi,coefi, "*")
	norm=sqrt(sum(H.theta*prod.coef))
	theta.normal[i,]=coefi/norm
}
Bspline.t0.theta<-splineDesign(delta.theta, t0, order.Bspline)
theta.t0<-theta.normal%*%t(Bspline.t0.theta)
pos<-which(theta.t0>0) 
THETA.seq.normalizado<-theta.normal[pos,]
num.norm <- nrow(THETA.seq.normalizado)
k.vec <- rep(0,num.norm)                                                                        
CV.app <- 0

if(is.null(n.core)) n.core<-availableCores(omit=1)
cl <- makeCluster(n.core)
registerDoParallel(cl)
m<-NULL
results <- foreach(m = 1:num.norm,.export="symsolve") %dopar% {
	cv.kseq <- 0  
	theta<-THETA.seq.normalizado[m,]
	theta.rec<-crossprod(t(Bspline.theta),theta) 
	Dmat.theta<-crossprod(Bspline,theta.rec)
	Theta.coef<-symsolve(Cmat,Dmat.theta)
	theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef})
	coef <-t(H.x%*%theta.x1);coef
	projec1<-rowSums(coef)
	semimetric <- outer(projec1,projec1,"-")
	norm.diff.0<-abs(semimetric)
	for(i in 1:n) {   
	    norm.order <- order(norm.diff.0[i,])
		zz <- sort(norm.diff.0[i,])[2:(kmax + 2)]
		bandwith <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
		z <- zz[ - (kmax + 1)]
		Zmat <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
		Umat <- Zmat/bandwith
		Kmat <- kernel(Umat)
		Kmat[col(Kmat) > row(Kmat)] <- 0
		ind.curves1 <- norm.order[2:(kmax + 1)]
		yind <- y[ind.curves1,]                 
		Ymat <- matrix(rep(yind, kmax), nrow = kmax, byrow = T)
		yhat1 <- rowSums(Ymat[k.seq,  ] * Kmat[k.seq,  ])/rowSums(Kmat[k.seq,  ])
		resid.kseq.2 <- (yhat1- y[i])^2
		cv.kseq <- cv.kseq + resid.kseq.2 		
	} 
	cv.kseq <- cv.kseq/(n*length.curve.y)
	index <- which.min(cv.kseq)
	k.opt.m <- k.seq[index]
	k.vec[m] <- k.opt.m
	CV.app[m] <- cv.kseq[index]
	for(j in 1:n) {  
		norm.diff.0j<-norm.diff.0[j,]
		norm.order <- order(norm.diff.0j)
		ind.curves2 <- norm.order[2:(k.opt.m + 2)]
		h<- sum(abs(norm.diff.0j[ind.curves2[k.opt.m:(k.opt.m + 1)]]))/2
		res.kernel <- kernel(norm.diff.0j[ind.curves2[ - (k.opt.m + 1)]]/h)	
		sum.res.kernel <- sum(res.kernel)
		yhat.cv[j,] <-ifelse(sum.res.kernel > 0, sum(y[ind.curves2[ - (k.opt.m + 1)],] * res.kernel)/sum.res.kernel,y[ind.curves2[1],])
		}
	estimated.Y[[m]] <- yhat.cv
	return(list(CV.app = CV.app[m], k.vec = k.vec[m], estimated.Y = estimated.Y[[m]]))		
}
stopCluster(cl)
CV.app <- sapply(results, function(x) x$CV.app)
m.opt <- order(CV.app)[1]
k.opt <- results[[m.opt]]$k.vec
theta.opt <- THETA.seq.normalizado[m.opt,]
theta.rec<-crossprod(t(Bspline.theta),theta.opt) 
Dmat.theta<-crossprod(Bspline,theta.rec)
Theta.coef<-symsolve(Cmat,Dmat.theta)
theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef})
coef <-t(H.x%*%theta.x1)
projec1<-rowSums(coef)
semimetric <- outer(projec1,projec1,"-")
norm.diff.1<-abs(semimetric)
H<-matrix(0,nrow=n,ncol=n)
for(i in 1:n){
	norm.diff.1i<-norm.diff.1[i,]
	norm.order <- order(norm.diff.1i)
	ind.curves2 <- norm.order[1:(k.opt + 1)]
	h <- sum(abs(norm.diff.1i[ind.curves2[k.opt:(k.opt + 1)]]))/2
	res.kernel <- kernel(norm.diff.1i[ind.curves2[ - (k.opt + 1)]]/h)
	sum.res.kernel <- sum(res.kernel)
	for(j in 1:k.opt){
		cur<-ind.curves2[j]
		if(sum.res.kernel > 0) {
			H[i,cur] <-res.kernel[j]/sum.res.kernel
		}                                                
	else H[i,cur] <- 1
	}
}
CV.opt<-CV.app[m.opt]
fitted.values<-H%*%y
res<-y-drop(fitted.values)
df<-sum(diag(H))
sigma.res<-sum(res^2)/(n-df)
r.2<-1-sum(res^2)/sum((y-mean(y))^2)
call<-match.call()
out<-list(fitted.values=fitted.values,residuals=res,theta.est=theta.opt, k.opt=k.opt, r.squared=r.2,var.res=sigma.res, df=df, yhat.cv=results[[m.opt]]$estimated.Y,CV.opt=CV.opt,CV.values=CV.app, H=H,m.opt=m.opt,theta.seq.norm=THETA.seq.normalizado, k.seq=k.seq,call=call,y=y,x=x,n=n,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,order.Bspline=order.Bspline,nknot.theta=nknot.theta)
class(out)<- "fsim.kNN"
out
}


