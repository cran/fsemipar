fsim.kNN.test<- function(x,y,x.test,y.test=NULL,
theta,order.Bspline=3,nknot.theta=3, 
k=4, 
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
y.hat2 <- matrix(0,J, ncol(y)) 
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 										  
for(j in 1:J) {       
	norm.diff <- semimetric.projec(data1=x, data2=x.test[j,], theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot = nknot, nknot.theta=nknot.theta)
	norm.order <- order(norm.diff)
	ind.curves2 <- norm.order[1:(k + 1)]
	h <- sum(abs(norm.diff[ind.curves2[k:(k + 1)]]))/2
	res.kernel <- kernel(norm.diff[ind.curves2[ - (k + 1)]]/h)
	res.kernel.mat <- matrix(rep(res.kernel, length.curve.y), k, length.curve.y, byrow=FALSE)
	sum.res.kernel <- sum(res.kernel)
	if(sum.res.kernel > 0) {
		if (length.curve.y>1)  y.hat2[j,] <- apply(y[ind.curves2[ - (k + 1)],] * res.kernel.mat,2,sum)/sum.res.kernel
		else  y.hat2[j,] <- sum(y[ind.curves2[ - (k + 1)],] * res.kernel)/sum.res.kernel
	}
	else y.hat2[j,] <- y[ind.curves2[1],]
}        
if (is.null(y.test)) MSEP<-NULL
else MSEP <- mean((y.hat2 - y.test)^2)
list(y.estimated.test=y.hat2, MSE.test=MSEP)
}

