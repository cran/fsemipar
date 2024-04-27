print.PVS.kernel<-function(x,...){
cat("*** MFPLM fitted using PVS with kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Bandwidth (h): ")
cat(x$h.opt)
cat("\n-wn: ")
cat(x$w.opt)
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


summary.PVS.kernel<-function(object,...){
x<-object
cat("*** MFPLM fitted using PVS with kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Bandwidth (h): ")
cat(x$h.opt)
cat("\n-wn: ")
cat(x$w.opt)
cat("\n-Linear coefficients (beta): \n")
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
cat("\n")
}


predict.PVS.kernel<- function(object,newdata.x=NULL,newdata.z=NULL,y.test=NULL,option=NULL, ...)
{
x<-object
if(is.null(newdata.x)|is.null(newdata.z)){
	y <- fitted(x)
	out<-y
}
else{
	if(is.null(option)) option<-1
	x.test <- newdata.x
	z.test<- newdata.z
	pred.LR.n <- as.matrix(z.test)%*%x$beta.est
    y.new<-x$y - as.matrix(x$z)%*%x$beta.est		
    y.new2<-x$y[x$train.2]-as.matrix(x$z[x$train.2,])%*%x$beta.est
    kind.of.semimetric <- paste("semimetric.", x$semimetric, sep = "")		
	if (option==1) {
		pred.FSIM.n <- fnp.kernel.test(y=y.new2,x=x$x[x$train.2,], x.test=x.test,y.test=y.test, kind.of.semimetric=kind.of.semimetric, h=x$h.opt, 
						q=x$q,kind.of.kernel=x$kind.of.kernel, range.grid=x$range.grid,nknot=x$nknot)
		pred.n <- pred.LR.n + pred.FSIM.n$y.estimate.test
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
		aux2 <-fnp.kernel.fit(y=y.new,x=x$x, kind.of.kernel=x$kind.of.kernel,kind.of.semimetric=kind.of.semimetric,min.quantile.h=x$min.quantile.h,
				 max.quantile.h=x$max.quantile.h,start.order.deriv.o.pca=x$q, end.order.deriv.o.pca=x$q,
				 h.seq=x$h.seq,num.h=x$num.h, range.grid=x$range.grid, nknot=x$nknot)
		h.opt.2 <- aux2$h.opt
		pred.FSIM.n.2 <- fnp.kernel.test(y=y.new,x=x$x, x.test=x.test, y.test=y.test,kind.of.semimetric=kind.of.semimetric, 
							h=h.opt.2, kind.of.kernel=x$kind.of.kernel, range.grid=x$range.grid, q=x$q,
							nknot=x$nknot)
		pred.n.2 <- pred.LR.n + pred.FSIM.n.2$y.estimate.test
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
	if ((option==3) & (length(x$indexes.beta.nonnull) != 0)) {
		aux <- sfpl.kernel.fit(x=x$x, z=x$z[,x$indexes.beta.nonnull], y=x$y,
				lambda.seq=0, nfolds=x$nfolds, seed=x$seed, min.q.h=x$min.quantile.h, 
				max.q.h=x$max.quantile.h, vn=x$vn, semimetric=x$semimetric,
				q=x$q, nknot=x$nknot, h.seq=x$h.seq, num.h=x$num.h, range.grid=x$range.grid,
				kind.of.kernel=x$kind.of.kernel,criterion=x$criterion, penalty=x$penalty, max.iter=x$max.iter)
		if(is.null(y.test)){
			MSEP.3<-NULL
			out<-fitted(aux)
		}
		else{
			MSEP3a<- predict(aux,newdata.x=x.test,newdata.z=z.test[,x$indexes.beta.nonnull],y.test=y.test,option=1)
			MSEP3b<- predict(aux,newdata.x=x.test,newdata.z=z.test[,x$indexes.beta.nonnull],y.test=y.test,option=2)    
				out<-list(y3a=MSEP3a$y,y3b=MSEP3b$y, MSEP3a=MSEP3a$MSEP.1, MSEP3b=MSEP3b$MSEP.2)
		}		 
	}
}  
out
}


plot.PVS.kernel<-function(x,ind=1:10,size=15,col1=1,col2=2,col3=4,option=0,...)
{
a<-x$range.grid[1]
b<-x$range.grid[2]
long_data <- as.data.frame(t(x$z[ind, ]))
Wavelength = seq(a, b, length.out = nrow(long_data))
long_data <- cbind(Wavelength, long_data)
long_data <- gather(long_data, key = "Series", value = "Value", -Wavelength)
matplot_colors <- 1:10
impact_points <- Wavelength[x$indexes.beta.nonnull]
g1=ggplot(long_data, aes(x = Wavelength, y = Value, group = Series, color = Series)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = rep(matplot_colors, length.out = length(unique(long_data$Series)))) +
  geom_vline(xintercept = impact_points, color = 1, linetype = "dashed",linewidth=1) +  
  theme_bw() +
  labs(x = "t", y = "", title = expression(paste(zeta[i],"(t)"))) +
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size),legend.position = "none")  
fit<-fitted(x)
mod <- lm(x$y ~ fit)
res <- residuals(mod)
y<-x$y
x_df <- data.frame(fit, y)
g2<-ggplot(x_df, aes(x = fit, y = y)) +
  geom_point(colour = col1,shape=1, size=5) + 
  geom_smooth(method=lm,formula=y~x,colour =col2,  linewidth = 1.5) +
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
if(option==0){
grid.arrange(g1,g2, g3,ncol =3)}
if(option==1){
grid.arrange(g1,ncol=1)}
if(option==2){
grid.arrange(g2,g3,ncol=2)}
}


