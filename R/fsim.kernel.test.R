fsim.kernel.test<- function(x,y,x.test,y.test, 
theta=theta, nknot.theta=3, order.Bspline=3, 
h=0.5,  
kind.of.kernel="quad", range.grid=NULL, nknot=NULL)
{
if (!is.matrix(x))  stop("x must contain a  matrix")
if (!is.matrix(x.test))  x.test <- t(as.matrix(x.test))
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
p <- ncol(x)
J <- nrow(x.test)
length.curve.y<-ncol(y) 
y.hat2 <- matrix(0,J, length.curve.y) 
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 					
for(j in 1:J) {   
	norm.diff <- semimetric.projec(data1=x, data2=x.test[j,], theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot = nknot, nknot.theta=nknot.theta)	
	res.kernel <- kernel(norm.diff/h)
	res.kernel[res.kernel<0] <- 0
	res.kernel[res.kernel>1] <- 0
	res.kernel.mat <- matrix(rep(res.kernel, length.curve.y), n, length.curve.y, byrow=FALSE)
	sum.res.kernel <- sum(res.kernel)
	if(sum.res.kernel > 0) {
		if (length.curve.y>1)  y.hat2[j,] <- apply(y * res.kernel.mat,2,sum)/sum.res.kernel
		else  y.hat2[j,] <- sum(y * res.kernel)/sum.res.kernel
	}
	else y.hat2[j,] <- y[order(norm.diff)[1],]
}       
if (is.null(y.test)) MSEP<-NULL
else MSEP <- mean((y.hat2 - y.test)^2)
list(y.estimated.test = y.hat2, MSE.test=MSEP)
}



