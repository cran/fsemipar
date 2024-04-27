fsim.kNN.fit.optim<- function(x,y, order.Bspline=3, nknot.theta=3, gamma=NULL, knearest=NULL,
min.knn=2, max.knn=NULL, step=NULL,  
kind.of.kernel="quad",range.grid=NULL,nknot=NULL,threshold=5e-3)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
p <- ncol(x)
length.curve.y<-ncol(y)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
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
if (is.null(gamma)) gamma <- c(rep(1, order.Bspline+nknot.theta)) 
prod.coef <-outer(gamma,gamma, "*")
norm<-sqrt(sum(H.theta*prod.coef))
theta<-gamma/norm
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
if (is.null(max.knn)) max.knn <- n%/%2
if (is.null(knearest)){
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
k.seq <- knearest
kmax <- max(k.seq)       
num.band <- length(k.seq)  
fsim.optim.kNN<-function(theta,y,Bspline.theta,H.theta,Cmat, Bspline,coef.mat1,H.x, k.opt,kernel)
{			
	prod.coef <-outer(theta,theta, "*")
	norm<-sqrt(sum(H.theta*prod.coef))
	theta<-theta/norm					   	
	theta.rec<-crossprod(t(Bspline.theta),theta) 
	Dmat.theta<-crossprod(Bspline,theta.rec)
	Theta.coef<-symsolve(Cmat,Dmat.theta)
	theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef})
	coef <-t(H.x%*%theta.x1);coef
	projec1<-rowSums(coef)
	projec2<-rowSums(coef)
	semimetric <- outer(projec1,projec2,"-")
	norm.diff.0<-abs(semimetric)
	n <- length(y)
	y.hat2 <- matrix(0,n,1) 
	k<-k.opt
	for(j in 1:n) {  
		norm.diff.0j<-norm.diff.0[j,]
		norm.order <- order(norm.diff.0j)
		ind.curves2 <- norm.order[2:(k + 2)] 
		h<- sum(abs(norm.diff.0j[ind.curves2[k:(k + 1)]]))/2
		res.kernel <- kernel(norm.diff.0j[ind.curves2[ - (k + 1)]]/h)	
		res.kernel.mat <- matrix(res.kernel,k, 1, byrow=FALSE)
		sum.res.kernel <- sum(res.kernel)
		y.hat2[j,] <-ifelse(sum.res.kernel > 0, sum(y[ind.curves2[ - (k + 1)],] * res.kernel)/sum.res.kernel,y[ind.curves2[1],])
		}
mean((y-y.hat2)^2)														
}
m <- 1
err2 <- 100
d.err <- 100                                                                                      
while(d.err >= threshold){
	m <- m + 1
	theta.rec<-crossprod(t(Bspline.theta),theta) 
	Dmat.theta<-crossprod(Bspline,theta.rec)
	Theta.coef<-symsolve(Cmat,Dmat.theta)
	theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef})
	coef <-t(H.x%*%theta.x1);coef
	projec1<-rowSums(coef)
	projec2<-rowSums(coef)
	semimetric <- outer(projec1,projec2,"-")
	norm.diff.0<-abs(semimetric)
	cv.kseq <- 0
	for (i in 1:n){
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
	k.opt<- k.seq[which.min(cv.kseq)]	
	theta.ini<-optim(par=theta,fn=fsim.optim.kNN, y=y,Bspline=Bspline,Bspline.theta=Bspline.theta,H.theta=H.theta,Cmat=Cmat,coef.mat1=coef.mat1,H.x=H.x,k.opt=k.opt,kernel=kernel)
	prod.coef <-outer(as.numeric(theta.ini$par),as.numeric(theta.ini$par), "*")
	norm<-sqrt(sum(H.theta*prod.coef))
	theta<-theta.ini$par/norm
	theta.rec<-crossprod(t(Bspline.theta),theta) 
	Dmat.theta<-crossprod(Bspline,theta.rec)
	Theta.coef<-symsolve(Cmat,Dmat.theta)
	theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef})
	coef <-t(H.x%*%theta.x1);coef
	projec1<-rowSums(coef)
	projec2<-rowSums(coef)
	semimetric <- outer(projec1,projec2,"-")
	norm.diff.0<-abs(semimetric)
	y.hat2 <- matrix(0,n,1) 
	for(j in 1:n) {  
		norm.diff.0j<-norm.diff.0[j,]
		norm.order <- order(norm.diff.0j)
		ind.curves2 <- norm.order[2:(k.opt + 2)] 
		h<- sum(abs(norm.diff.0j[ind.curves2[k.opt:(k.opt + 1)]]))/2
		res.kernel <- kernel(norm.diff.0j[ind.curves2[ - (k.opt + 1)]]/h)	
		res.kernel.mat <- matrix(res.kernel,k.opt, 1, byrow=FALSE)
		sum.res.kernel <- sum(res.kernel)
		y.hat2[j,] <-ifelse(sum.res.kernel > 0, sum(y[ind.curves2[ - (k.opt + 1)],] * res.kernel)/sum.res.kernel,y[ind.curves2[1],])
		}
	CV.tot<-mean((y-y.hat2)^2)	
	err2[m]<-CV.tot/var(y)
    if(err2[m] - err2[m-1] < 0) d.err <- abs(err2[m]-err2[m-1]) else d.err <- threshold + 1      
    if(m >= 20) d.err <- threshold*0.9
}
H<-matrix(0,nrow=n,ncol=n)
for(i in 1:n){
	norm.diff.1i<-norm.diff.0[i,]
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
df<-sum(diag(H))
fitted.values<-H%*%y
res<-y-drop(fitted.values)
sigma.res<-sum(res^2)/(n-df)
r.2<-1-sum(res^2)/sum((y-mean(y))^2)
call<-match.call()
out<-list(fitted.values=fitted.values,residuals=res,theta.est=theta, k.opt=k.opt, r.squared=r.2,var.res=sigma.res, df=df, err2=err2[-1], H=H, k.seq=k.seq, CV.opt=CV.tot,CV.kseq=cv.kseq,
		call=call,y=y,x=x,n=n,
		kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
		order.Bspline=order.Bspline,nknot.theta=nknot.theta)
class(out)<- "fsim.kNN"
out
}


