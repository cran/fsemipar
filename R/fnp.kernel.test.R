fnp.kernel.test <- function(x,y, x.test,y.test=NULL,
kind.of.semimetric="semimetric.deriv", q=2, ext.inf=1, ext.sup=3, 
h=0.5, 
kind.of.kernel="quad", range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(x.test))  x.test<- t(as.matrix(x.test))
kernel <- get(kind.of.kernel)
semimetric <- get(kind.of.semimetric)
n <- nrow(x)
p <- ncol(x)
J <- nrow(x.test)
if (!is.matrix(y))  y <- as.matrix(y)
if (is.null(range.grid)) range.grid <- c(1,p)
length.curve.y<-ncol(y)
y.hat2 <- matrix(0,J, length.curve.y) 
nknot.m <- nknot
if (is.null(nknot)) {
	if (kind.of.semimetric=="semimetric.interv") nknot.m <- (p - 0 - 3 - 1)%/%2 
	else  if (kind.of.semimetric=="semimetric.deriv") nknot.m <- (p - q - 3 - 1)%/%2 
}
for(j in 1:J) {
	if (p==1) norm.diff <- abs(x - x.test[j,])                                           
	else if (kind.of.semimetric=="semimetric.deriv") norm.diff <- semimetric(data1=x, data2=x.test[j,], q=q, range.grid=range.grid, nknot = nknot.m)
	else if (kind.of.semimetric=="semimetric.interv") norm.diff <- semimetric(data1=x, data2=x.test[j,], interv=c(ext.inf,ext.sup), range.grid=range.grid, nknot = nknot.m)                     
	else if (kind.of.semimetric=="semimetric.pca")  norm.diff <- semimetric(data1=x, data2=x.test[j,], q=q)                     
	res.kernel <- kernel(norm.diff/h)
	res.kernel[res.kernel<0] <- 0
	res.kernel[res.kernel>1] <- 0
	res.kernel.mat <- matrix(rep(res.kernel, length.curve.y), n, length.curve.y, byrow=FALSE)
	sum.res.kernel <- sum(res.kernel)
	if(sum.res.kernel > 0) {
		if (length.curve.y>1)  y.hat2[j,] <- colSums(y * res.kernel.mat)/sum.res.kernel
		else  y.hat2[j,] <- sum(y * res.kernel)/sum.res.kernel
	}
	else y.hat2[j,] <- y[order(norm.diff)[1],]
} # for (j       
if(!is.null(y.test))MSEP <- mean((y.hat2 - y.test)^2)
else MSEP<-NULL
estimated.Y <- y.hat2     
list(y.estimate.test = y.hat2, MSE.test=MSEP,y.test=y.test)
} 



