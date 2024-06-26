fsim.kernel.fit.fixedtheta<- function(y,x, 
norm.diff,min.quantile.h=0.05, max.quantile.h=0.5, h.seq=NULL, num.h=10,kind.of.kernel)
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
norm.diff<-as.matrix(norm.diff)
cv.hseq<-numeric(num.h)
if (is.null(h.seq)) {
    Semimetric.0 <- norm.diff[row(norm.diff) > col(norm.diff)]
    h.seq <- quantile(Semimetric.0, seq(min.quantile.h, max.quantile.h, length = num.h))
}
for (j in 1:num.h) {
		h <- h.seq[j]
		res.kernel<-kernel(norm.diff/h)
        res.kernel[res.kernel<0] <- 0
	    res.kernel[res.kernel>1] <- 0	
        sum.res.kernel <-colSums(res.kernel)
        yhat1<-res.kernel%*%y/sum.res.kernel
		input.num<-y[apply(norm.diff,2,order)[2,]]
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
res.kernel<-kernel(norm.diff/h.opt)
res.kernel[res.kernel<0] <- 0
res.kernel[res.kernel>1] <- 0	
sum.res.kernel <-colSums(res.kernel)
yhat.cv<-res.kernel%*%y/sum.res.kernel
input.num<-y[apply(norm.diff,2,order)[2,]]
yhat.cv[is.na(yhat.cv)]<-input.num[is.na(yhat.cv)]
h.min.opt.max<-c(min(h.seq),h.opt,max(h.seq))
estimated.Y<-yhat.cv	
list(h.opt=h.opt, h.seq=h.seq, h.min.opt.max=h.min.opt.max, yhat.cv = estimated.Y, y=y, CV=CV.app)
} 



