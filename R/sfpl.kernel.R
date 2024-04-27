print.sfpl.kernel<-function(x,...){
cat("*** SFPL fitted using penalized least squares combined with kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Bandwidth (h): ")
cat(x$h.opt)
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est)
cat("\n-Number of non-zero linear coefficients: ")
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


summary.sfpl.kernel<-function(object,...){
x<-object
cat("*** SFPL fitted using penalized least squares combined with kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Bandwidth (h): ")
cat(x$h.opt)
cat("\n-Linear coefficients (beta): ")
cat(x$beta.est)
cat("\n-Number of non-zero linear coefficients: ")
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


predict.sfpl.kernel<- function(object,newdata.x=NULL,newdata.z=NULL,y.test=NULL,option=NULL, ...)
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
		kind.of.semimetric <- paste("semimetric.", object$semimetric, sep = "")
		pred.FNP.n <- fnp.kernel.test(y=y.new,x=object$x, y.test=y.test, x.test=x.test, kind.of.semimetric=kind.of.semimetric, q=object$q, h=object$h.opt, 
						kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, nknot=object$nknot)
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
		kind.of.semimetric <- paste("semimetric.", object$semimetric, sep = "")
		pred.FNP.n <- fnp.kernel.fit.test(y=y.new, x=object$x,y.test=y.test, x.test=x.test, kind.of.semimetric=kind.of.semimetric, kind.of.kernel=object$kind.of.kernel, 
						min.quantile.h=object$min.quantile.h, max.quantile.h= object$max.quantile.h, h.seq=object$h.seq, num.h=object$num.h, start.order.deriv.o.pca=object$q, 
						end.order.deriv.o.pca=object$q, range.grid=object$range.grid, nknot=object$nknot)
		pred.n <- pred.LR.n + pred.FNP.n$y.estimate.test
		y2<-pred.n
		h.option2<-pred.FNP.n$h.opt
		if(is.null(y.test)){
			MSEP.2<-NULL
			out<-list(y=y2,h.opt.2=h.option2)
		}
		else{
			MSEP.2 <-mean((y2 - y.test)^2)
			out<-list(y=y2,MSEP.2=MSEP.2,h.opt.2=h.option2)
		} 		 
    }
}  
out
}


plot.sfpl.kernel<-function(x,size=15,col1=1,col2=2,col3=4,...)
{
fit<-fitted(x)
mod <- lm(x$y ~ fit)
res <- residuals(mod)
y<-x$y
x_df <- data.frame(fit, y)
g2<-ggplot(x_df, aes(x = fit, y = y)) +
  geom_point(colour = col1,shape=1, size=5) + 
  geom_smooth(method=lm,formula=y~x,colour =col2, linewidth= 1.5) +
  labs(x = "Fitted values", y = "y", title = "Response vs Fitted values") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
x_df2 <- data.frame(fit, res)
g3<-ggplot(x_df2, aes(x = fit, y = res)) +
  geom_point(colour = col1,shape=1,size=5) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = 1, linewidth = 1) +
  geom_smooth(method=loess,formula=y~x,linewidth=1.5,col=col3)+
  labs(x = "Fitted values", y = "Residuals", title = "Residuals vs Fitted Values") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
grid.arrange(g2, g3,ncol = 2)
}

