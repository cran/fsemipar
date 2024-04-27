print.PVS.kNN<-function(x,...){
cat("*** MFPLM fitted using PVS with kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Number of neighbours (k): ")
cat(x$k.opt)
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


summary.PVS.kNN<-function(object,...){
x<-object
cat("*** MFPLM fitted using PVS with kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Number of neighbours (k): ")
cat(x$k.opt)
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


predict.PVS.kNN<- function(object,newdata.x=NULL,newdata.z=NULL,y.test=NULL,option=NULL,knearest.n=object$knearest, min.knn.n=object$min.knn, max.knn.n=object$max.knn.n, step.n=object$step, ...)
{
if(is.null(newdata.x)|is.null(newdata.z)){
	y <- fitted(object)
	out<-y}
else{
	if(is.null(option)) option<-1
	x.test <- newdata.x
	z.test<- newdata.z
	pred.LR.n <- as.matrix(z.test)%*%object$beta.est
	y.new<-object$y - as.matrix(object$z)%*%object$beta.est		
	y.new2<-object$y[object$train.2]-as.matrix(object$z[object$train.2,])%*%object$beta.est
	kind.of.semimetric <- paste("semimetric.", object$semimetric, sep = "")
	if (option==1) {
		pred.FNP.n <- fnp.kNN.test(y=y.new2,x=object$x[object$train.2,], x.test=x.test,y.test=y.test, k=object$k.opt, kind.of.semimetric=kind.of.semimetric,q=object$q,
			kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid,nknot=object$nknot)
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
		aux2 <-fnp.kNN.fit(y=y.new,x=object$x, kind.of.semimetric=kind.of.semimetric, kind.of.kernel=object$kind.of.kernel,min.knn=min.knn.n,max.knn=max.knn.n,
				knearest=knearest.n,step=step.n, range.grid=object$range.grid, nknot=object$nknot,  start.order.deriv.o.pca=object$q, end.order.deriv.o.pca=object$q)
		k.opt.2 <- aux2$k.opt
		pred.FSIM.n.2 <- fnp.kNN.test(y=y.new,x=object$x, x.test=x.test, y.test=y.test, k=k.opt.2,kind.of.semimetric=kind.of.semimetric,q=object$q,
							kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, nknot=object$nknot)
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
	if ((option==3) & (length(object$indexes.beta.nonnull) != 0)) {
		aux <- sfpl.kNN.fit(x=object$x, z=object$z[,object$indexes.beta.nonnull], y=object$y,nknot=object$nknot, semimetric=object$semimetric,q=object$q,
				lambda.seq=0,vn=object$vn, nfolds=object$nfolds, seed=object$seed, knearest=knearest.n, min.knn=min.knn.n, max.knn=max.knn.n, step=step.n, 
				range.grid=object$range.grid,kind.of.kernel=object$kind.of.kernel, criterion=object$criterion, penalty=object$penalty,
				max.iter=object$max.iter)
		if(is.null(y.test)){
			MSEP.3<-NULL
			out<-fitted(aux)
		}
		else{
			MSEP3a<- predict(aux,newdata.x=x.test,newdata.z=z.test[,object$indexes.beta.nonnull],y.test=y.test,option=1)
			MSEP3b<- predict(aux,newdata.x=x.test,newdata.z=z.test[,object$indexes.beta.nonnull],y.test=y.test,option=2)    
			out<-list(y3a=MSEP3a$y,y3b=MSEP3b$y, MSEP3a=MSEP3a$MSEP.1,MSEP3b=MSEP3b$MSEP.2)      
		}		 
	}
    if (option==4) {
		pred.FNP <- fnp.kNN.fit.test.loc(x=object$x,y=y.new, x.test=x.test,y.test=y.test, kind.of.semimetric=kind.of.semimetric, knearest=knearest.n, min.knn=min.knn.n, 
						max.knn=max.knn.n, step=step.n, kind.of.kernel=object$kind.of.kernel, start.order.deriv.o.pca=object$q, end.order.deriv.o.pca=object$q, min.leng.interv=NULL, 
						max.leng.interv=NULL, range.grid=object$range.grid, nknot=object$nknot)
		pred <- pred.LR.n + pred.FNP$y.estimate.test
		y4<-pred
		if(is.null(y.test)){
			MSEP.4<-NULL
			out<-y4
		}
		else{
			MSEP.4 <-mean((y4 - y.test)^2)
			out<-list(y=y4,MSEP.4=MSEP.4)
		}
	}            
}  
out
}


plot.PVS.kNN<-function(x,ind=1:10,size=15,col1=1,col2=2,col3=4,option=0,...)
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

