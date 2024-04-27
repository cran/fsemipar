print.fsim.kernel<-function(x,...){
cat("*** FSIM fitted using kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Bandwidth (h): ")
cat(x$h.opt)
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


summary.fsim.kernel<-function(object,...){
x<-object
cat("*** FSIM fitted using kernel estimation with Nadaraya-Watson weights ***\n")
cat("\n-Call: ")
print(x$call)
cat("\n-Bandwidth (h): ")
cat(x$h.opt)
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


predict.fsim.kernel<- function(object,newdata=NULL,y.test=NULL, ...)
{
if(is.null(newdata)){
	y <- fitted(object)
	out<-y
}
else{
	x.test <- newdata
	res<- fsim.kernel.test(y=object$y,x=object$x, x.test=x.test,y.test=y.test, theta=object$theta.est, h=object$h.opt, 
	kind.of.kernel=object$kind.of.kernel, range.grid=object$range.grid, order.Bspline=object$order.Bspline, 
	nknot=object$nknot, nknot.theta=object$nknot.theta)
	y<-res$y.estimated.test
	if(is.null(y.test)) out<-y
	else out<-list(y=y,MSEP=res$MSE.test)
}
out
}


plot.fsim.kernel<-function(x, size=15,col1=1, col2=2,...)
{
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
theta_df <- data.frame(x.t, theta.rec)
g1<-ggplot(theta_df, aes(x = x.t, y = theta.rec)) +
  geom_line(linewidth = 1.5, color = col1) +
  labs(x = "range.grid X", y = "", title = expression(widehat(theta)[0])) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))  
x.hat.theta <- projec(data = x$x, theta = x$theta.est, range.grid = x$range.grid, order.Bspline = x$order.Bspline, nknot = x$nknot, nknot.theta = x$nknot.theta)
y.hat <- fitted(x)
vec <- cbind(x.hat.theta, y.hat)
vec2 <- vec[order(x.hat.theta),]
regression_df <- data.frame(x.hat.theta, x$y)
g2<-ggplot(regression_df, aes(x = x.hat.theta, y = x$y)) +
  geom_point(color = "black",size=5,shape=1) +
  geom_line(data = data.frame(vec2), aes(x = vec[,1], y = vec[,2]), color = col2, linewidth = 1.5) +
  labs(x = expression(paste("<", widehat(theta)[0], ",", "X", ">")), y = "", title =expression(paste(widehat(r),"(","<", widehat(theta)[0], ",", "X", ">",")")) ) +
  theme_bw()+
theme(plot.title = element_text(hjust = 0.5,size=size),axis.title.x=element_text(size=size))

grid.arrange(g1, g2, ncol = 2)
}