print.FASSMR.kNN<-function(x,...){
cat("*** MFPLSIM fitted using FASSMR combined with kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Number of neighbours (k): ")
cat(x$k.opt)
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


summary.FASSMR.kNN<-function(object,...){
cat("*** MFPLSIM fitted using FASSMR combined with kNN estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(object$call)
cat("\n-Bandwidth (h): ")
cat(object$k.opt)
cat("\n-wn: ")
cat(object$w.opt)
cat("\n-Theta coefficients in the B-spline basis: ")
cat(object$theta.est)
cat("\n-Linear coefficients (beta): \n")
cat(object$beta.est)
cat("\n-Number of non-zero linear coefficients: ")
cat(length(object$indexes.beta.nonnull))
cat("\n-Indexes non-zero beta-coefficients: ")
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


predict.FASSMR.kNN<- function(object,newdata.x=NULL,newdata.z=NULL,y.test=NULL,option=NULL, ...)
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
		pred.FSIM.n <- fsim.kNN.test(y=y.new,x=object$x, x.test=x.test,y.test=y.test, theta=object$theta.est, k=object$k.opt, 
						kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, order.Bspline=object$order.Bspline, 
						nknot=object$nknot, nknot.theta=object$nknot.theta)
		pred.n <- pred.LR.n + pred.FSIM.n$y.estimated.test
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
		aux2 <-fsim.kNN.fit.fixedtheta(y=y.new,x=object$x,norm.diff=object$norm.diff, min.knn=object$min.knn,max.knn=object$max.knn,knearest=object$knearest,
				step=object$step,kind.of.kernel=object$kind.of.kernel)
		k.opt.2 <- aux2$k.opt
		pred.FSIM.n.2 <- fsim.kNN.test(y=y.new,x=object$x, x.test=x.test, y.test=y.test, theta=object$theta.est, k=k.opt.2, kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, 
							order.Bspline=object$order.Bspline, nknot=object$nknot, nknot.theta=object$nknot.theta)
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


plot.FASSMR.kNN<-function(x,ind=1:10,size=15,col1=1,col2=2,col3=4,option=0,...)
{
a<-x$range.grid[1]
b<-x$range.grid[2]
long_data <- as.data.frame(t(x$z[ind, ]))
Wavelength = seq(a, b, length.out = nrow(long_data))
long_data <- cbind(Wavelength, long_data)
long_data <- gather(long_data, key = "Series", value = "Value", -Wavelength)
matplot_colors <- 1:10
impact_points <- Wavelength[x$indexes.beta.nonnull]
g1=ggplot(long_data, aes(x = Wavelength, y = Value, Series, color = Series)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = rep(matplot_colors, length.out = length(unique(long_data$Series)))) +
  geom_vline(xintercept = impact_points, color = 1, linetype = "dashed",linewidth=1) +  
  theme_bw() +
  labs(x = "t", y = "", title = expression(paste(zeta[i],"(t)"))) +
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size),legend.position = "none")  
  
THETA<-x$theta.est
nknot.theta<-x$nknot.theta
order.Bspline<-x$order.Bspline
x.t <- seq(a, b, length=ncol(x$x))
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
Bspline.theta<-splineDesign(delta.theta,x.t,order.Bspline)
theta.rec<-Bspline.theta%*%THETA 
theta_df <- data.frame(x.t, theta.rec)
g2<-ggplot(theta_df, aes(x = x.t, y = theta.rec)) +
  geom_line(linewidth = 1.5, color = col1) +
  labs(x = "range.grid X", y = "", title = expression(widehat(theta)[0])) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
fit<-fitted(x)
mod <- lm(x$y ~ fit)
res <- residuals(mod)
y<-x$y
x_df <- data.frame(fit, y)
g3<-ggplot(x_df, aes(x = fit, y = y)) +
  geom_point(colour = col1,shape=1, size=5) + 
  geom_smooth(method=lm,formula=y~x,colour =col2,  linewidth = 1.5) +
  labs(x = "Fitted values", y = "y", title = "Response vs Fitted values") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
x_df2 <- data.frame(fit, res)
g4<-ggplot(x_df2, aes(x = fit, y = res)) +
  geom_point(colour = col1,shape=1,size=5) + 
  geom_hline(yintercept = 0, linetype = "dashed", colour = 1, linewidth = 1) +
  geom_smooth(method=loess,formula=y~x,linewidth=1.5,col=col3)+
  labs(x = "Fitted Values", y = "Residuals", title = "Residuals vs Fitted Values") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))
if(option==0){
grid.arrange(g1,g2,g3,g4,ncol =2,nrow=2)}
if(option==1){
grid.arrange(g1,ncol=1)}
if(option==2){
grid.arrange(g2,ncol=1)}
if(option==3){
grid.arrange(g3,g4,ncol=2)}
}





