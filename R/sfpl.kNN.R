
print.sfpl.kNN<-function(x,...){
cat("*** SFPL model fitted using penalized least squares combined with kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Number of neighbours (k): ")
cat(x$k.opt)
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est)
cat("\n-Number of non-zero beta-coefficients: ")
cat(length(x$indexes.beta.nonnull))
cat("\n-Indexes non-zero beta-coefficients: ")
cat(x$indexes.beta.nonnull)
cat("\n-Lambda: ")
cat(x$lambda.opt)
cat("\n-IC: ")
cat(x$IC)
cat("\n-Penalty: ")
cat(x$penalty)
cat("\n-Criterion: ")
cat(x$criterion)
cat("\n-vn: ")
cat(x$vn.opt)
cat("\n")
}


summary.sfpl.kNN<-function(object,...){
x<-object
cat("*** SFPL model fitted using penalized least squares combined with kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Number of neighbours (k): ")
cat(x$k.opt)
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est)
cat("\n-Number of non-zero beta-coefficients: ")
cat(length(x$indexes.beta.nonnull))
cat("\n-Indexes non-zero beta-coefficients: ")
cat(x$indexes.beta.nonnull)
cat("\n-Lambda: ")
cat(x$lambda.opt)
cat("\n-IC: ")
cat(x$IC)
cat("\n-Penalty: ")
cat(x$penalty)
cat("\n-Criterion: ")
cat(x$criterion)
cat("\n-vn: ")
cat(x$vn.opt)
cat("\n")
}


predict.sfpl.kNN<- function(object,newdata.x=NULL,newdata.z=NULL,y.test=NULL,option=NULL, ...)
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
	resY.T.n <- object$y - as.matrix(object$z)%*%object$beta.est 																												  
    if (option==1) {
		kind.of.semimetric <- paste("semimetric.", object$semimetric, sep = "")
		pred.FNP.n <- fnp.kNN.test(y=resY.T.n,x=object$x,x.test=x.test,y.test=y.test,kind.of.semimetric=kind.of.semimetric, q=object$q,
						k=object$k.opt, kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, nknot=object$nknot)
		pred.n <- pred.LR.n + pred.FNP.n$y.estimate.test
		y1<-pred.n
		if(is.null(y.test)){
			MSEP.1<-NULL
			out<-y1
		}
		else {
			MSEP.1 <-  mean((y1 - y.test)^2)
			out<-list(y=y1,MSEP.1=MSEP.1)  			
		}
	}
	if (option==2) {
		pred.FNP.n <- fnp.kNN.fit.test(y=resY.T.n, x=object$x, x.test=x.test, y.test=y.test, semimetric=object$semimetric, knearest=object$knearest,
						kind.of.kernel=object$kind.of.kernel, start.order.deriv.o.pca=object$q, end.order.deriv.o.pca=object$q, min.leng.interv=NULL,
						max.leng.interv=NULL, range.grid=object$range.grid, nknot=object$nknot)
		pred.n <- pred.LR.n + pred.FNP.n$y.estimate.test
		y2<-pred.n
		knn.option2<- pred.FNP.n$k.opt
		if(is.null(y.test)){
			MSEP.2<-NULL
			out<-y2
		}
		else {
			MSEP.2 <-  mean((y2 - y.test)^2)
			out<-list(y=y2,MSEP.2=MSEP.2,k.opt.2=knn.option2)  	
		}	
	}
	if (option==3) { 
		kind.of.semimetric <- paste("semimetric.", object$semimetric, sep = "")
		pred.FNP.n <- fnp.kNN.fit.test.loc(y=resY.T.n,x=object$x, y.test=y.test,x.test=x.test, kind.of.semimetric=kind.of.semimetric, knearest=object$knearest,
						kind.of.kernel=object$kind.of.kernel, start.order.deriv.o.pca=object$q, end.order.deriv.o.pca=object$q, min.leng.interv=NULL, 
						max.leng.interv=NULL, range.grid=object$range.grid, nknot=object$nknot)
		pred.n <- pred.LR.n + pred.FNP.n$y.estimate.test
		y3<-pred.n
		if(is.null(y.test)){
			MSEP.3<-NULL
			out<-y3
		}
		else{
			MSEP.3 <-  mean((y3 - y.test)^2)
			out<-list(y=y3,MSEP.3=MSEP.3) 
		} 
	}			
}             
out
}


plot.sfpl.kNN<-function(x,cex.axis=1.5,cex.lab=1.5,cex=2,col=1,cex.main=1.5,...)
{
oldpar <- par(no.readonly = TRUE)    
on.exit(par(oldpar))
par(mfrow=c(1,2))
plot(x$fitted.values,x$y,xlab="Fitted values", ylab="y",cex.lab=cex.lab,cex.axis=cex.axis,cex=cex,col=col,cex.main=cex.main,main="Response vs fitted values")
mod<-lm(x$y~x$fitted.values)
abline(mod, col=2,lwd=2)
plot(x$fitted.values,x$residuals,xlab="Fitted values", ylab="Residuals",cex.lab=cex.lab,cex.axis=cex.axis,cex=cex,cex.main=cex.main,col=col,main="Residuals vs fitted values")
abline(h=0,lty=2,lwd=2)
}
