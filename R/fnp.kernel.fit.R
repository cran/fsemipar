fnp.kernel.fit<- function(x,y, 
kind.of.semimetric="semimetric.deriv",start.order.deriv.o.pca=NULL, end.order.deriv.o.pca=NULL, min.leng.interv=NULL, max.leng.interv=NULL,  
min.quantile.h=0.05, max.quantile.h=0.5, h.seq=NULL, num.h=10,
kind.of.kernel="quad",range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if (!is.matrix(x))  stop("x must contain a matrix")
kernel <- get(kind.of.kernel)
semimetric <- get(kind.of.semimetric)
n <- nrow(x)
if (!(is.null(h.seq))) num.h <- length(h.seq)
h.seq.aux <- list()
p <- ncol(x)
if (!is.matrix(y))  y <- as.matrix(y)
length.curve.y<-ncol(y)
estimated.Y <- list()
yhat.cv <- matrix(0,n,length.curve.y)
if (p==1) { 
	start.order <- 1
	end.order <- 1
}
else 
	if ((kind.of.semimetric=="semimetric.deriv") | (kind.of.semimetric=="semimetric.pca")){		       
	if ( (kind.of.semimetric=="semimetric.deriv") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 0
	if ( (kind.of.semimetric=="semimetric.deriv") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- 2
	if ((kind.of.semimetric=="semimetric.pca") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 1
	if ((kind.of.semimetric=="semimetric.pca") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- p
	start.order <- start.order.deriv.o.pca
	end.order <- end.order.deriv.o.pca
}
else
	if (kind.of.semimetric=="semimetric.interv"){
	if (is.null(min.leng.interv)) min.leng.interv <- 1
	if (is.null(max.leng.interv)) max.leng.interv <- p-1
	start.order <- 1
	end.order <- (max.leng.interv-min.leng.interv +1)*(2*p-max.leng.interv-min.leng.interv)/2
}
else stop("kind.of.semimetric must be equal to: semimetric.deriv, semimetric.interv or semimetric.pca")
if (is.null(range.grid)) range.grid <- c(1,p)
nknot.m <- nknot
if (is.null(nknot)) {
	if (kind.of.semimetric=="semimetric.interv") nknot.m <- (p - 0 - 3 - 1)%/%2 
}
num.norm <- end.order - start.order +1 	  
CV.app <- 0
h.seq.app <- list()
h.vec <- rep(0,num.norm)									
for(m in start.order:end.order) {
	if (is.null(nknot)) {
	if (kind.of.semimetric=="semimetric.deriv") nknot.m <- (p - m - 3 - 1)%/%2 
	}
	a <- p-min.leng.interv
	if (kind.of.semimetric=="semimetric.interv") {
		for (i in 1:(max.leng.interv-min.leng.interv+1)) 
		if ( ((i-1)*(2*a+2-i)/2 + 1 <= m) & (m <= (i*(2*a+1-i)/2)) ) {		                  
			ext.inf <- m - (i-1)*(2*a+2-i)/2
			ext.sup <- ext.inf + min.leng.interv + i -1
		} 
	}		
	h.seq.m <- h.seq
	num.h.m <- num.h
	if (is.null(h.seq)) {
		if (p==1) norm.diff.0 <- abs(outer(x, x, "-"))			
		else if (kind.of.semimetric=="semimetric.deriv") norm.diff.0 <- semimetric(data1=x, data2=x, q=m, range.grid=range.grid, nknot = nknot.m)
		else if (kind.of.semimetric=="semimetric.interv") norm.diff.0 <- semimetric(data1=x, data2=x, interv=c(ext.inf,ext.sup), range.grid=range.grid, nknot = nknot.m)					  
		else if (kind.of.semimetric=="semimetric.pca")  norm.diff.0 <- semimetric(data1=x, data2=x, q=m)
		norm.diff.00 <- norm.diff.0[row(norm.diff.0) > col(norm.diff.0)]
		h.seq.m <- quantile(norm.diff.00, seq(min.quantile.h, max.quantile.h, length = num.h))
		h.seq.m <- h.seq.m[h.seq.m>0]
		num.h.m <- length(h.seq.m)
	}
	cv.hseq <- rep(0,num.h.m)
	h.seq.app[[m + 1 - start.order]] <- h.seq.m
	for(i in 1:n) {       			   
		for (j in 1:num.h.m) {
			h <- h.seq.m[j]
			yhat1 <- fnp.kernel.test(x=x[-i,],y=y[-i],y.test=y[i], x.test=x[i,], kind.of.semimetric=kind.of.semimetric, q=m, h=h, kind.of.kernel=kind.of.kernel, range.grid=range.grid, nknot=nknot)$y.estimate.test
			if (length.curve.y==1) resid.hseq.2 <- (yhat1- y[i])^2
			else resid.hseq.2 <- sum((yhat1- y[i,])^2)
			cv.hseq[j] <- cv.hseq[j] + resid.hseq.2
			} # for (j		                 
	} # for (i
	cv.hseq <- cv.hseq/n
	index <- order(cv.hseq)[1]
	h.opt.m <- h.seq.m[index]
	h.vec[m + 1 - start.order] <- h.opt.m
	CV.app[m + 1 - start.order] <- cv.hseq[index]
	for (i in 1:n) 
		yhat.cv[i,] <- fnp.kernel.test(x=x[-i,], y=y[-i], x.test=x[i,], y.test=y[i], kind.of.semimetric=kind.of.semimetric, q=m, h=h.opt.m, kind.of.kernel=kind.of.kernel, range.grid=range.grid, nknot=nknot)$y.estimate.test
	estimated.Y[[m + 1 - start.order]] <- yhat.cv	 
} 
m.opt <- order(CV.app)[1] +  start.order - 1
h.opt <- h.vec[m.opt + 1 - start.order]
h.seq.opt <- h.seq.app[[m.opt + 1 - start.order]]
if (kind.of.semimetric=="semimetric.interv") {		
	a <- p-min.leng.interv   	
	for (i in 1:(max.leng.interv-min.leng.interv+1)) 
	if ( ((i-1)*(2*a+2-i)/2 + 1 <= m.opt) & (m.opt <= (i*(2*a+1-i)/2)) ) {
		ext.inf <- m.opt - (i-1)*(2*a+2-i)/2
		ext.sup <- ext.inf + min.leng.interv + i -1
	}
	list(m.opt=m.opt, h.opt=h.opt, h.seq.opt=h.seq.opt, interv.opt=c(ext.inf, ext.sup),  y.estimate.app = estimated.Y[[order(CV.app)[1]]], y=y,
	CV=CV.app)
} 
	else	
	list(m.opt=m.opt, h.opt=h.opt, h.seq.opt=h.seq.opt, y.estimate.app = estimated.Y[[order(CV.app)[1]]], y=y,
	CV=CV.app)
} 



