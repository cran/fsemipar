print.FASSMR.kernel<-function(x,...){
cat("*** MFPLSIM fitted using FASSMR combined with kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Bandwidth (h): ")
cat(x$h.opt)
cat("\n-wn: ")
cat(x$w.opt)
cat("\n-Theta coefficients in the B-spline basis: ")
cat(x$theta.est)
cat("\n-Linear coefficients (beta): \n")
cat(x$beta.est)
cat("\n-Number of non-zero beta-coefficients: ")
cat(length(x$indexes.beta.nonnull))
cat("\n-Indexes non-zero linear coefficients: ")
cat(x$indexes.beta.nonnull)
cat("\n-Lambda: ")
cat(x$lambda.opt)
cat("\n-IC: ")
cat(x$IC)
cat("\n-Penalty: ")
cat(x$penalty)
cat("\n-Criterion: ")
cat(x$criterion)
cat("\n")
}


summary.FASSMR.kernel<-function(object,...){
cat("*** MFPLSIM fitted using FASSMR combined with kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(object$call)
cat("\n-Bandwidth (h): ")
cat(object$h.opt)
cat("\n-wn: ")
cat(object$w.opt)
cat("\n-Theta coefficients in the B-spline basis: ")
cat(object$theta.est)
cat("\n-Linear coefficients (beta): \n")
cat(object$beta.est)
cat("\n-Number of non-zero linear coefficients: ")
cat(length(object$indexes.beta.nonnull))
cat("\n-Indexes of non-zero beta-coefficients: ")
cat(object$indexes.beta.nonnull)
cat("\n-Lambda: ")
cat(object$lambda.opt)
cat("\n-IC: ")
cat(object$IC)
cat("\n-Penalty: ")
cat(object$penalty)
cat("\n-Criterion: ")
cat(object$criterion)
cat("\n")
}


predict.FASSMR.kernel<- function(object,newdata.x=NULL,newdata.z=NULL,y.test=NULL,option=NULL, ...)
{
if(is.null(newdata.x)|is.null(newdata.z)){
	y <- fitted(object)
	out<-y
}
else{
	if(is.null(option)) option<-1
	x.test <- newdata.x
	z.test<- newdata.z
	pred.LR.n <- as.matrix(z.test)%*%object$beta.est
	y.new<-object$y - as.matrix(object$z)%*%object$beta.est		
	if (option==1) {
		pred.FSIM.n <- fsim.kernel.test(y=y.new,x=object$x, x.test=x.test,y.test=y.test, theta=object$theta.est, h=object$h.opt, 
						kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, order.Bspline=object$order.Bspline, 
						nknot=object$nknot, nknot.theta=object$nknot.theta)
		pred.n <- pred.LR.n + pred.FSIM.n$y.estimated.test
		y1<-pred.n
		if(is.null(y.test)){
			MSEP.1<-NULL
			out<-y1
		}
		else{
			MSEP.1 <-  mean((y1 - y.test)^2)
			out<-list(y=y1,MSEP.1=MSEP.1)       
			}
	}
	if (option==2) {
		aux2 <-fsim.kernel.fit.fixedtheta(y=y.new,x=object$x, theta=object$theta.est, kind.of.kernel=object$kind.of.kernel,min.quantile.h=object$min.quantile.h,max.quantile.h=object$max.quantile.h,
				h.seq=object$h.seq,num.h=object$num.h, range.grid=object$range.grid, order.Bspline=object$order.Bspline, nknot=object$nknot, nknot.theta=object$nknot.theta)
		h.opt.2 <- aux2$h.opt
		pred.FSIM.n.2 <- fsim.kernel.test(y=y.new,x=object$x, x.test=x.test, y.test=y.test, theta=object$theta.est, h=h.opt.2, kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, 
							order.Bspline=object$order.Bspline,nknot=object$nknot, nknot.theta=object$nknot.theta)
		pred.n.2 <- pred.LR.n + pred.FSIM.n.2$y.estimated.test
		y2<-pred.n.2
		if(is.null(y.test)){
			MSEP.2<-NULL
			out<-y2
		}
		else{
			MSEP.2 <-mean((y2 - y.test)^2)
			out<-list(y=y2,MSEP.2=MSEP.2)
		} 		 
	}
}  
out
}


plot.FASSMR.kernel<-function(x,cex.axis=1.5,cex.lab=1.5,cex=2,col=1,cex.main=1.5,...)
{
oldpar <- par(no.readonly = TRUE)    
on.exit(par(oldpar))
par(mfrow=c(1,3))
plot(x$fitted.values,x$y,xlab="Fitted values", ylab="y",cex.lab=cex.lab,cex.axis=cex.axis,cex=cex,col=col,cex.main=cex.main,main="Response vs fitted values")
mod<-lm(x$y~x$fitted.values)
abline(mod, col=2,lwd=2)
plot(x$fitted.values,x$residuals,xlab="Fitted values", ylab="Residuals",cex.lab=cex.lab,cex.axis=cex.axis,cex=cex,cex.main=cex.main,col=col,main="Residuals vs fitted values")
abline(h=0,lty=2)
THETA<-x$theta.est
a<-x$range.grid[1]
b<-x$range.grid[2]
nknot.theta<-x$nknot.theta
order.Bspline<-x$order.Bspline
nm<-length(x$x[1,])
x.t <- seq(a, b, length=nm)
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
Bspline.theta<-splineDesign(delta.theta,x.t,order.Bspline)
theta.rec<-Bspline.theta%*%THETA 
plot(x.t,theta.rec,type="l",xlim=c(a,b),ylab="", xlab="range.grid X",cex.main=cex.main,main=expression(widehat(theta)[0]),lwd=2,cex.lab=cex.lab,cex.axis=cex.axis,col=col)
}



