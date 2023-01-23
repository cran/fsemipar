fsim.kNN.fit.fixedtheta<- function(y,x, 
theta, order.Bspline=3,nknot.theta=3,
min.knn=2, max.knn=NULL, knearest=NULL, step=NULL, 
kind.of.kernel="quad",range.grid=NULL,nknot=NULL)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
p <- ncol(x)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(max.knn)) max.knn <- n%/%2
if (is.null(knearest)) {
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
num.knn <- length(knearest)
estimated.Y <- list()
length.curve.y<-ncol(y)
yhat.cv <- matrix(0,n,length.curve.y)
dim.base.theta <- order.Bspline + nknot.theta
cv.kseq <- rep(0,num.knn)
for(i in 1:n) {       
	for (j in 1:num.knn) {
		k <- knearest[j]
		yhat1 <- fsim.kNN.test(y=y[-i],x=x[-i,], y.test=y[i],x.test=x[i,], theta=theta, k=k, kind.of.kernel=kind.of.kernel, range.grid=range.grid, 
							   order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)$y.estimated.test
		if (length.curve.y==1) resid.kseq.2 <- (yhat1- y[i])^2
		else resid.kseq.2 <- sum((yhat1- y[i,])^2)
		cv.kseq[j] <- cv.kseq[j] + resid.kseq.2
	} 	
} 
cv.kseq <- cv.kseq/(n*length.curve.y)
index <- order(cv.kseq)[1]
k.opt <- knearest[index]
knn.min.opt.max <- c(min(knearest), k.opt, max(knearest))
CV.app <- cv.kseq[index]
for (i in 1:n) 
	yhat.cv[i,] <- fsim.kNN.test(y=y[-i],x=x[-i,], y.test=y[i],x.test=x[i,], theta=theta, k=k.opt, kind.of.kernel=kind.of.kernel, range.grid=range.grid,
								 order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)$y.estimated.test
estimated.Y <- yhat.cv	 
list(k.opt=k.opt, knearest=knearest, knn.min.opt.max=knn.min.opt.max, yhat.cv= estimated.Y, y=y, CV=CV.app)
} 



