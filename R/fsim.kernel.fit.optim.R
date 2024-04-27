fsim.kernel.fit.optim<- function(x,y,
nknot.theta=3,order.Bspline=3, gamma=NULL,
min.q.h=0.05, max.q.h=0.5, h.seq=NULL, num.h=10, 
kind.of.kernel="quad", range.grid=NULL, nknot=NULL,threshold=5e-3)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
min.quantile.h<-min.q.h
max.quantile.h<-max.q.h
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
if(is.null(gamma)) gamma <- c(rep(1, order.Bspline+nknot.theta)) 
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

fsim.optim<-function(theta,y,Bspline.theta,H.theta,Cmat,Bspline,coef.mat1,H.x, h.opt,kernel)
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
	res.kernel<-kernel(norm.diff.0/h.opt)
	res.kernel[res.kernel<0] <- 0
	res.kernel[res.kernel>1] <- 0	
	sum.res.kernel <-colSums(res.kernel)
	yhat2<-res.kernel%*%y/sum.res.kernel
	input.num<-y[apply(norm.diff.0,2,order)[2,]]
	yhat2[is.na(yhat2)]<-input.num[is.na(yhat2)]		     
	den<-1-kernel(0)/sum.res.kernel
	dif<-((y-yhat2)/den)^2	
	dif[is.na(dif)]<-(y[is.na(dif)]-input.num[is.na(dif)])^2
	sum(dif)														
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
	norm.diff.00 <- norm.diff.0[row(norm.diff.0) > col(norm.diff.0)]
	h.seq.m <- quantile(norm.diff.00, seq(min.quantile.h, max.quantile.h, length = num.h))
	h.seq.m <- h.seq.m[h.seq.m>0]
	num.h.m <- length(h.seq.m)
	cv.hseq <- rep(0,num.h.m) 
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
	h.opt<- h.seq.m[which.min(cv.hseq)]
    theta.ini<-optim(par=theta,fn=fsim.optim, y=y,Bspline=Bspline,Bspline.theta=Bspline.theta,H.theta=H.theta,Cmat=Cmat,coef.mat1=coef.mat1,H.x=H.x,h.opt=h.opt,kernel=kernel)
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
	res.kernel<-kernel(norm.diff.0/h.opt)
    res.kernel[res.kernel<0] <- 0
	res.kernel[res.kernel>1] <- 0	
    sum.res.kernel <-colSums(res.kernel)
    yhat1<-res.kernel%*%y/sum.res.kernel
	input.num<-y[apply(norm.diff.0,2,order)[2,]]
    yhat1[is.na(yhat1)]<-input.num[is.na(yhat1)]		     
	den<-1-kernel(0)/sum.res.kernel
	res.cv<-((y-yhat1)/den)^2	
    res.cv[is.na(res.cv)]<-(y[is.na(res.cv)]-input.num[is.na(res.cv)])^2
	CV.tot<-mean(res.cv)		 
    err2[m] <- CV.tot/var(y)
    if(err2[m] - err2[m-1] < 0) d.err <- abs(err2[m]-err2[m-1]) else d.err <- threshold + 1      
    if(m >= 20) d.err <- threshold*0.9
}
H<-res.kernel/sum.res.kernel
H[is.na(H)]<-1
df<-sum(diag(H))
fitted.values<-H%*%y
res<-y-drop(fitted.values)
sigma.res<-sum(res^2)/(n-df)
r.2<-1-sum(res^2)/sum((y-mean(y))^2)
call<-match.call()
out<-list(fitted.values=fitted.values, residuals=res,theta.est=theta, h.opt=h.opt,r.squared=r.2,var.res=sigma.res,df=df,
		  err=err2[-1],H=H, h.seq=h.seq.m, CV.opt=CV.tot,CV.hseq=cv.hseq,
		  call=call,y=y,x=x,n=n,
		  kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
		  order.Bspline=order.Bspline,nknot.theta=nknot.theta)
class(out)<-"fsim.kernel"
out
} 

