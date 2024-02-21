fsim.kernel.fit.fixedtheta<- function(y,x, 
theta, order.Bspline=3, nknot.theta=3,
min.quantile.h=0.05, max.quantile.h=0.5, h.seq=NULL, num.h=10,
kind.of.kernel="quad",range.grid,nknot)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
p <- ncol(x)
if (!(is.null(h.seq))) num.h <- length(h.seq)
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
for (j in 1:num.h) {
		h <- h.seq[j]
		res.kernel<-kernel(norm.diff.0/h)
        res.kernel[res.kernel<0] <- 0
	    res.kernel[res.kernel>1] <- 0	
        sum.res.kernel <-colSums(res.kernel)
        yhat1<-res.kernel%*%y/sum.res.kernel
		input.num<-y[apply(norm.diff.0,2,order)[2,]]
        yhat1[is.na(yhat1)]<-input.num[is.na(yhat1)]		     
		den<-1-kernel(0)/sum.res.kernel
	    dif<-((y-yhat1)/den)^2	
        dif[is.na(dif)]<-(y[is.na(dif)]-input.num[is.na(dif)])^2
	    cv.hseq[j] <-sum(dif)		 
} 	
index <- which.min(cv.hseq)
cv.hseq <- cv.hseq/n
h.opt<- h.seq[index]
CV.app<- cv.hseq[index]
res.kernel<-kernel(norm.diff.0/h.opt)
res.kernel[res.kernel<0] <- 0
res.kernel[res.kernel>1] <- 0	
sum.res.kernel <-colSums(res.kernel)
yhat.cv<-res.kernel%*%y/sum.res.kernel
input.num<-y[apply(norm.diff.0,2,order)[2,]]
yhat.cv[is.na(yhat.cv)]<-input.num[is.na(yhat.cv)]
h.min.opt.max<-c(min(h.seq),h.opt,max(h.seq))
estimated.Y<-yhat.cv	
list(h.opt=h.opt, h.seq=h.seq, h.min.opt.max=h.min.opt.max, yhat.cv = estimated.Y, y=y, CV=CV.app)
} 



