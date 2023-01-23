
print.fsim.kNN<-function(x,...){
cat("*** FSIM fitted using kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Number of neighbours (k): ")
cat(x$k.opt)
cat("\n-Theta coefficients in the B-spline basis: ")
cat(x$theta.est)
cat("\n-CV: ")
cat(x$CV.opt)
cat("\n-R squared: ")
cat(x$r.squared)
cat("\n-Residual variance: ")
cat(x$var.res)
cat(" on ")
cat(x$n-x$df)
cat(" degrees of freedom")
cat("\n")
}


summary.fsim.kNN<-function(object,...){
x<-object
cat("*** FSIM fitted using kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Number of neighbours (k): ")
cat(x$k.opt)
cat("\n-Theta coefficients in the B-spline basis: ")
cat(x$theta.est)
cat("\n-CV: ")
cat(x$CV.opt)
cat("\n-R squared: ")
cat(x$r.squared)
cat("\n-Residual variance: ")
cat(x$var.res)
cat(" on ")
cat(x$n-x$df)
cat(" degrees of freedom")
cat("\n")
}


predict.fsim.kNN<- function(object,newdata=NULL,y.test=NULL, ...)
{
if(is.null(newdata)){
	y <- fitted(object)
	out<-y
}
else{
	x.test <- newdata
	res<- fsim.kNN.test(y=object$y,x=object$x, x.test=x.test,y.test=y.test, theta=object$theta.est, k=object$k.opt, 
			kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, order.Bspline=object$order.Bspline, 
			nknot=object$nknot, nknot.theta=object$nknot.theta)
	y<-res$y.estimated.test
	if(is.null(y.test)) out<-y
	else out<-list(y=y,MSEP=res$MSE.test)
}
out
}


plot.fsim.kNN<-function(x,cex.axis=1.5,cex.lab=1.5,cex=2,col=1,cex.main=1.5,...)
{
oldpar <- par(no.readonly = TRUE)    
on.exit(par(oldpar))
par(mfrow=c(1,2))
THETA<-x$theta.est
a<-x$range.grid[1]
b<-x$range.grid[2]
nknot.theta<-x$nknot.theta
order.Bspline<-x$order.Bspline
x.t <- seq(a, b, length=ncol(x$x))
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
Bspline.theta<-splineDesign(delta.theta,x.t,order.Bspline)
theta.rec<-Bspline.theta%*%THETA 
plot(x.t,theta.rec,type="l",xlim=c(a,b),ylab="", xlab="range.grid X",main=expression(widehat(theta)[0]),cex.axis=cex.axis, cex.lab=cex.lab, cex=cex,col=col, cex.main=cex.main)
x.hat.theta=projec(data=x$x, theta=x$theta.est, range.grid=x$range.grid,order.Bspline=x$order.Bspline, nknot=x$nknot,nknot.theta=x$nknot.theta)
y.hat=fitted(x)
vec=cbind(x.hat.theta,y.hat)
vec2=vec[order(x.hat.theta),] 
plot(x.hat.theta,x$y,xlab=expression(paste("<",widehat(theta)[0],",","X",">")),ylab="",main="Regression fit",cex.axis=cex.axis, cex.lab=cex.lab, cex=cex,col=col, cex.main=cex.main)
lines(vec2[,1],vec2[,2],type="l",col=2,lwd=2)
}


