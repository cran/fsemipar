fsim.kernel.fit<- function(x,y,
seed.coeff=c(-1,0,1), order.Bspline=3, nknot.theta=3,
min.q.h=0.05, max.q.h=0.5, h.seq=NULL, num.h=10, 
kind.of.kernel="quad", range.grid=NULL, nknot=NULL,n.core=NULL)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
min.quantile.h<-min.q.h
max.quantile.h<-max.q.h
if (!(is.null(h.seq))) num.h <- length(h.seq)
p <- ncol(x)
length.curve.y<-ncol(y)
estimated.Y <- list()
yhat.cv <- matrix(0,n,length.curve.y)
if (is.null(range.grid)) range.grid <- c(1,p)
t0 <- mean(range.grid)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
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
theta.rec<-crossprod(t(Bspline.theta),t(THETA.seq.normalizado)) 
Dmat.theta<-crossprod(Bspline,theta.rec)
Theta.coef<-symsolve(Cmat,Dmat.theta)	
num.norm <- nrow(THETA.seq.normalizado)
h.seq.app <- list()
h.vec <- rep(0,num.norm)									
CV.app <- 0  
if(is.null(n.core)) n.core<-availableCores(omit=1)
cl <- makeCluster(n.core)
registerDoParallel(cl)
m<-NULL
results <- foreach(m = 1:num.norm) %dopar% {  
	h.seq.m <- h.seq
	num.h.m <- num.h
	theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef[,m]})
	coef <-t(H.x%*%theta.x1)
	projec1<-rowSums(coef)
	semimetric <- outer(projec1,projec1,"-")
	norm.diff.0<-abs(semimetric)
	if (is.null(h.seq)){ 
		norm.diff.00 <- norm.diff.0[row(norm.diff.0) > col(norm.diff.0)]
		h.seq.m <- quantile(norm.diff.00, seq(min.quantile.h, max.quantile.h, length = num.h))
		h.seq.m <- h.seq.m[h.seq.m>0]
		num.h.m <- length(h.seq.m)
		cv.hseq <- rep(0,num.h.m) 
	}   
	h.seq.app[[m]] <- h.seq.m
	for (j in 1:num.h.m) {
		h <- h.seq.m[j]
		res.kernel<-kernel(norm.diff.0/h)
        res.kernel[res.kernel<0] <- 0
	    res.kernel[res.kernel>1] <- 0	
        sum.res.kernel <-colSums(res.kernel)
        yhat1<-res.kernel%*%y/sum.res.kernel
		input.num<-y[apply(norm.diff.0,2,order)[2,]]
        yhat1[is.na(yhat1)]<-input.num[is.na(yhat1)]		     
		den<-1-kernel(0)/sum.res.kernel
	    dif<-((y-yhat1)/den)^2	
        dif[is.na(dif)]<-(y[is.na(dif)]-input.num[is.na(dif)])^2
	    cv.hseq[j] <-sum(dif)		 
	} 	
	cv.hseq <- cv.hseq/n
	index <- which.min(cv.hseq)
	h.opt.m <- h.seq.m[index]
	h.vec[m] <- h.opt.m
	CV.app[m] <- cv.hseq[index]
	res.kernel<-kernel(norm.diff.0/h.opt.m)
	res.kernel[res.kernel<0] <- 0
	res.kernel[res.kernel>1] <- 0	
	sum.res.kernel <-colSums(res.kernel)
	yhat.cv<-res.kernel%*%y/sum.res.kernel
	input.num<-y[apply(norm.diff.0,2,order)[2,]]
	yhat.cv[is.na(yhat.cv)]<-input.num[is.na(yhat.cv)]
	estimated.Y[[m]] <- yhat.cv
	return(list(CV.app = CV.app[m], h.vec = h.vec[m], h.seq.app = h.seq.app[[m]], estimated.Y = estimated.Y[[m]]))
} 
stopCluster(cl)
CV.app <- sapply(results, function(x) x$CV.app)
m.opt <- order(CV.app)[1]
h.opt <- results[[m.opt]]$h.vec
h.seq.opt <- results[[m.opt]]$h.seq.app
theta.opt <- THETA.seq.normalizado[m.opt,]
theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef[,m.opt]})
coef <-t(H.x%*%theta.x1)
projec1<-rowSums(coef)
semimetric <- outer(projec1,projec1,"-")
norm.diff.1<-abs(semimetric)
H<-matrix(0,nrow=n,ncol=n)
res.kernel <- kernel(norm.diff.1/h.opt)
res.kernel[res.kernel<0] <- 0
res.kernel[res.kernel>1] <- 0
sum.res.kernel <- colSums(res.kernel)
H<-res.kernel/sum.res.kernel
CV.opt<-CV.app[m.opt]
fitted.values<-H%*%y
res<-y-drop(fitted.values)
df<-sum(diag(H))
sigma.res<-sum(res^2)/(n-df)
r.2<-1-sum(res^2)/sum((y-mean(y))^2)
call<-match.call()
out<-list(fitted.values=fitted.values, residuals=res,theta.est=theta.opt, h.opt=h.opt,r.squared=r.2,var.res=sigma.res,df=df,yhat.cv=results[[m.opt]]$estimated.Y,CV.opt=CV.opt,CV.values=CV.app,H=H,m.opt=m.opt, theta.seq.norm=THETA.seq.normalizado,  h.seq=h.seq.opt, 
call=call,y=y,x=x,n=n,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
order.Bspline=order.Bspline,nknot.theta=nknot.theta)
class(out)<-"fsim.kernel"
out
} 



