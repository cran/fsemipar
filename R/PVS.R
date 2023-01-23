print.PVS<-function(x,...){
cat("*** MLM fitted using PVS algorithm ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Intercept (beta_0): ")
cat(x$beta.est[1])
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est[-1])
cat("\n-Number of non-zero beta-coefficients: ")
cat(length(x$indexes.beta.nonnull[[1]]))
cat("\n-Indexes non-zero beta-coefficients: ")
cat(x$indexes.beta.nonnull[[1]])
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


summary.PVS<-function(object,...){
x<-object
cat("*** MLM fitted using PVS algorithm ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Intercept (beta_0): ")
cat(x$beta.est[1])
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est[-1])
cat("\n-Number of non-zero beta-coefficients: ")
cat(length(x$indexes.beta.nonnull[[1]]))
cat("\n-Indexes non-zero beta-coefficients: ")
cat(x$indexes.beta.nonnull[[1]])
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


predict.PVS<- function(object,newdata=NULL,y.test=NULL, ...)
{
if(is.null(newdata)){
	y <- fitted(object)
	out<-y
}
else{
	z.test <- newdata
	onez.test<- cbind(rep(1, nrow(z.test)), z.test)
	pred.LR <- as.matrix(onez.test)%*%object$beta.est
	MSEP<-  mean((pred.LR - y.test)^2)		
	if(is.null(y.test)) out<-pred.LR
	else out<-list(y=pred.LR,MSEP=MSEP)
}
out
}


plot.PVS<-function(x,cex.axis=1.5,cex.lab=1.5,cex=2,col=1,cex.main=1.5,...)
{
oldpar <- par(no.readonly = TRUE)    
on.exit(par(oldpar))
par(mfrow=c(1,2))
plot(x$fitted.values,x$y,xlab="Fitted values", ylab="y",cex.lab=cex.lab,cex.axis=cex.axis,cex=cex,col=col,cex.main=cex.main,main="Response vs fitted values")
mod<-lm(x$y~x$fitted.values)
abline(mod, col=2,lwd=2)
plot(x$fitted.values,x$residuals,xlab="Fitted values", ylab="Residuals",cex.lab=cex.lab,cex.axis=cex.axis,cex=cex,cex.main=cex.main,col=col,main="Residuals vs fitted values")
abline(h=0,lty=2)
}
