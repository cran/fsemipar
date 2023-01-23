fsim.kernel.fit.fixedtheta<- function(y,x, 
theta, order.Bspline=3, nknot.theta=3,
min.quantile.h=0.05, max.quantile.h=0.5, h.seq=NULL, num.h=10,
kind.of.kernel="quad",range.grid=NULL,nknot=NULL)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
p <- ncol(x)
if (!(is.null(h.seq))) num.h <- length(h.seq)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
estimated.Y <- list()
length.curve.y<-ncol(y)
yhat.cv <- matrix(0,n,length.curve.y)
dim.base.theta <- order.Bspline + nknot.theta
if (is.null(h.seq)) { 
	norm.diff.0 <- semimetric.projec(data1=x, data2=x, theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)				
	norm.diff.00 <- norm.diff.0[row(norm.diff.0) > col(norm.diff.0)]
	h.seq <- quantile(norm.diff.00, seq(min.quantile.h, max.quantile.h, length = num.h))
	h.seq <- h.seq[h.seq>0]
	num.h <- length(h.seq)
	cv.hseq <- rep(0,num.h) 
}
for(i in 1:n) {              
	for (j in 1:num.h) {
		h <- h.seq[j]
		yhat1 <- fsim.kernel.test(y=y[-i],x=x[-i,], y.test=y[i],x.test=x[i,], theta=theta, h=h, kind.of.kernel=kind.of.kernel, range.grid=range.grid, 
								  order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)$y.estimated.test
		if (length.curve.y==1) resid.hseq.2 <- (yhat1- y[i])^2
		else resid.hseq.2 <- sum((yhat1- y[i,])^2)
		cv.hseq[j] <- cv.hseq[j] + resid.hseq.2
	} 		                 
} 
cv.hseq <- cv.hseq/(n*length.curve.y)
index <- order(cv.hseq)[1]
h.opt <- h.seq[index]
h.min.opt.max <- c(min(h.seq), h.opt, max(h.seq))
CV.app <- cv.hseq[index]
for (i in 1:n) 
	yhat.cv[i,] <- fsim.kernel.test(y=y[-i],x=x[-i,], y.test=y[i],x.test=x[i,], theta=theta, h=h.opt,  kind.of.kernel=kind.of.kernel, range.grid=range.grid, 
									order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)$y.estimated.test
estimated.Y <- yhat.cv 	
list(h.opt=h.opt, h.seq=h.seq, h.min.opt.max=h.min.opt.max, yhat.cv = estimated.Y, y=y, CV=CV.app)
} 



