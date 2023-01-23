fnp.kNN.test<- function(y,x, y.test=NULL,x.test, 
kind.of.semimetric="semimetric.deriv", q=2,  ext.inf=1, ext.sup=3, 
k=4,
kind.of.kernel="quad", range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if (!is.matrix(x))  stop("x must be a matrix")
if (!is.matrix(x.test))  x.test<- t(as.matrix(x.test))
if (is.null(kind.of.kernel)) kind.of.kernel <- "quad"

kernel <- get(kind.of.kernel)
semimetric <- get(kind.of.semimetric)
n <- nrow(x)
p <- ncol(x)
J <- nrow(x.test)

if (is.null(range.grid)) range.grid <- c(1,p)
if (!is.matrix(y))  y <- as.matrix(y)

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
estimated.Y <- y.hat2     
list(y.estimate.test = y.hat2, MSE.test=MSEP)
} 



