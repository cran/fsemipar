fsim.kNN.fit.fixedtheta.com<- function(y,x, 
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
kmax<-max(knearest)
k.seq <- knearest
num.knn <- length(knearest)
estimated.Y <- list()
length.curve.y<-ncol(y)
yhat.cv <- matrix(0,n,length.curve.y)
dim.base.theta <- order.Bspline + nknot.theta
cv.kseq <- rep(0,num.knn)
norm.diff.0 <- semimetric.projec(data1=x, data2=x, theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)
for(i in 1:n) {   
	    norm.order <- order(norm.diff.0[i,])
		zz <- sort(norm.diff.0[i,])[2:(kmax + 2)]
		bandwith <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
		z <- zz[ - (kmax + 1)]
		Zmat <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
		Umat <- Zmat/bandwith
		Kmat <- kernel(Umat)
		Kmat[col(Kmat) > row(Kmat)] <- 0
		ind.curves1 <- norm.order[2:(kmax + 1)]
		yind <- y[ind.curves1,]                 
		Ymat <- matrix(rep(yind, kmax), nrow = kmax, byrow = T)
		yhat1 <- rowSums(Ymat[k.seq,  ] * Kmat[k.seq,  ])/rowSums(Kmat[k.seq,  ])
		resid.kseq.2 <- (yhat1- y[i])^2
		cv.kseq <- cv.kseq + resid.kseq.2 		
	} 
	cv.kseq <- cv.kseq/(n*length.curve.y)
	index <- which.min(cv.kseq)
	k.opt <- knearest[index]
	CV.app <- cv.kseq[index]
	knn.min.opt.max <- c(min(knearest), k.opt, max(knearest))
	for(j in 1:n) {  
		norm.diff.0j<-norm.diff.0[j,]
		norm.order <- order(norm.diff.0j)
		ind.curves2 <- norm.order[2:(k.opt + 2)]
		h<- sum(abs(norm.diff.0j[ind.curves2[k.opt:(k.opt+ 1)]]))/2
		res.kernel <- kernel(norm.diff.0j[ind.curves2[ - (k.opt + 1)]]/h)	
		sum.res.kernel <- sum(res.kernel)
		yhat.cv[j,] <-ifelse(sum.res.kernel > 0, sum(y[ind.curves2[ - (k.opt + 1)],] * res.kernel)/sum.res.kernel,y[ind.curves2[1],])
		}
estimated.Y <- yhat.cv 
list(k.opt=k.opt, knearest=knearest, knn.min.opt.max=knn.min.opt.max, yhat.cv= estimated.Y, y=y, CV=CV.app)
} 


