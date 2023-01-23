fnp.kNN.fit<- function(x=x,y=y, 
kind.of.semimetric="semimetric.deriv",  start.order.deriv.o.pca=NULL, end.order.deriv.o.pca=NULL, min.leng.interv=NULL, max.leng.interv=NULL,
knearest=NULL, min.knn=2, max.knn=NULL, step=NULL,  
kind.of.kernel="quad",range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if (is.null(kind.of.kernel)) kind.of.kernel <- "quad"
kernel <- get(kind.of.kernel)
if (is.null(kind.of.semimetric)) kind.of.semimetric <- "semimetric.deriv"
semimetric <- get(kind.of.semimetric)
if (!is.matrix(x)) x<-as.matrix(x)
if (!is.matrix(y)) y<-as.matrix(y)
n <- nrow(x)
p <- ncol(x)
length.curve.y<-ncol(y)
if (is.null(max.knn)) max.knn <- (n %/% 2)
estimated.Y <- list()
yhat.cv <- matrix(0,n,length.curve.y) 
if (p==1) { 
	start.order <- 1
	end.order <- 1
}
else if ((kind.of.semimetric=="semimetric.deriv") | (kind.of.semimetric=="semimetric.pca")){		       
	if ( (kind.of.semimetric=="semimetric.deriv") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 0
	if ( (kind.of.semimetric=="semimetric.deriv") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- 2
	if ((kind.of.semimetric=="semimetric.pca") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 1
	if ((kind.of.semimetric=="semimetric.pca") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- p
	start.order <- start.order.deriv.o.pca
	end.order <- end.order.deriv.o.pca
}
else if (kind.of.semimetric=="semimetric.interv"){
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
if (is.null(min.knn)) min.knn <- 2
if (is.null(max.knn)) max.knn <- n%/%2
if (is.null(knearest)) {
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
k.seq <- knearest
kmax <- max(k.seq)           
CV.app <- 0
MSEP <- 0
num.band <- length(k.seq)       
k.vec <- rep(0,num.norm)									       
for(m in start.order:end.order) {
	cv.kseq <- 0	
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
	for(i in 1:n) {                      
		if (p==1) norm.diff <- abs(x- x[i,])            
		else if (kind.of.semimetric=="semimetric.deriv") norm.diff <- semimetric(data1=x, data2=x[i,], q=m, range.grid=range.grid, nknot = nknot.m)
		else if (kind.of.semimetric=="semimetric.interv") norm.diff <- semimetric(data1=x, data2=x[i,], interv=c(ext.inf,ext.sup), range.grid=range.grid, nknot = nknot.m)                      
		else if (kind.of.semimetric=="semimetric.pca")  norm.diff <- semimetric(data1=x, data2=x[i,], q=m)
		norm.order <- order(norm.diff)
		zz <- sort(norm.diff)[2:(kmax + 2)]
		bandwith <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
		z <- zz[ - (kmax + 1)]
		Zmat <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
		Umat <- Zmat/bandwith
		Kmat <- kernel(Umat)
		Kmat[col(Kmat) > row(Kmat)] <- 0
		ind.curves1 <- norm.order[2:(kmax + 1)]
		yind <- y[ind.curves1,]
		if (length.curve.y==1)  {                       
			Ymat <- matrix(rep(yind, kmax), nrow = kmax, byrow = T)
			yhat1 <- apply(Ymat[k.seq,  ] * Kmat[k.seq,  ], 1, sum)/apply(Kmat[k.seq,  ], 1, sum)
			resid.kseq.2 <- (yhat1- y[i])^2
			cv.kseq <- cv.kseq + resid.kseq.2  
		}
		else {
			Karr <- array(t(Kmat), c(dim(t(Kmat)),length.curve.y))
			Karr <- aperm(Karr, c(1,3,2))
			Yarr <- array(yind, c(dim(yind),kmax))
			yhat1 <- (t(apply(Yarr[,,k.seq] * Karr[,,k.seq], 2:3, sum)))/matrix(apply(Kmat[k.seq,], 1, sum),length(k.seq),length.curve.y)
			resid.kseq.2 <- apply((yhat1- matrix(rep(y[i,],nrow(yhat1)), dim(yhat1), byrow=TRUE))^2, 1, sum)             
			cv.kseq <- cv.kseq + resid.kseq.2
		}                  
	} 
	cv.kseq <- cv.kseq/(n*length.curve.y)
	index <- order(cv.kseq)[1]
	k.opt.m <- k.seq[index]
	k.vec[m + 1 - start.order] <- k.opt.m
	CV.app[m + 1 - start.order] <- cv.kseq[index]
	for (i in 1:n) 
		yhat.cv[i,] <- fnp.kNN.test(x=x[-i,],y=y[-i],x.test=x[i,],y.test=y[i], kind.of.semimetric=kind.of.semimetric, q=m, k=k.opt.m, kind.of.kernel=kind.of.kernel, range.grid=range.grid, nknot=nknot)$y.estimate.test
	estimated.Y[[m + 1 - start.order]] <- yhat.cv	 
} 
m.opt <- order(CV.app)[1] +  start.order - 1
k.opt <- k.vec[m.opt + 1 - start.order]
if (kind.of.semimetric=="semimetric.interv") {		
	a <- p-min.leng.interv   	
	for (i in 1:(max.leng.interv-min.leng.interv+1)) 
	if ( ((i-1)*(2*a+2-i)/2 + 1 <= m.opt) & (m.opt <= (i*(2*a+1-i)/2)) ) {
		ext.inf <- m.opt - (i-1)*(2*a+2-i)/2
		ext.sup <- ext.inf + min.leng.interv + i -1
	}
	list(m.opt=m.opt, k.opt=k.opt, interv.opt=c(ext.inf, ext.sup),  Y.estimate.app = estimated.Y[[order(CV.app)[1]]], Y.app=y,
	CV.app=CV.app)
} 
	else	
	list(m.opt=m.opt, k.opt=k.opt,  y.estimate.app = estimated.Y[[order(CV.app)[1]]], y=y,
	CV.app=CV.app)
} 

