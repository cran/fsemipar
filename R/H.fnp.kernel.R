H.fnp.kernel<-function(x,h,
semimetric="deriv",q=2,ext.inf=1, ext.sup=3,
kind.of.kernel="quad",range.grid=NULL,nknot=NULL)
{
if (!is.matrix(x)) x<-as.matrix(x)
kind.of.semimetric <- paste("semimetric.", semimetric, sep = "")
p<-ncol(x)
if (is.null(range.grid)) range.grid <- c(1,p)
nknot.m <- nknot
if (is.null(nknot)) {
	if (kind.of.semimetric=="semimetric.interv") nknot.m <- (p - 0 - 3 - 1)%/%2 
	else  if (kind.of.semimetric=="semimetric.deriv") nknot.m <- (p - q - 3 - 1)%/%2 
}
n<-nrow(x)
kernel <- get(kind.of.kernel)
H<-matrix(0,nrow=n,ncol=n)
for(i in 1:n){
	if (p==1) norm.diff <- abs(x - x[i,])
	else if (kind.of.semimetric=="semimetric.deriv") norm.diff <- semimetric.deriv(data1=x[i,], data2=x, q=q, range.grid=range.grid, nknot = nknot.m)
	else if (kind.of.semimetric=="semimetric.interv") norm.diff <- semimetric.interv(data1=x[i,], data2=x, interv=c(ext.inf,ext.sup), range.grid=range.grid, nknot = nknot.m)
	else if (kind.of.semimetric=="semimetric.pca")  norm.diff <- semimetric.pca(data1=x[i,], data2=x, q=q)
	res.kernel <- kernel(norm.diff/h)
	res.kernel[res.kernel<0] <- 0
	res.kernel[res.kernel>1] <- 0
	sum.res.kernel <- sum(res.kernel)
	if(sum.res.kernel > 0) {
		H[i,] <-res.kernel/sum.res.kernel
	}                                                
	else H[i,order(norm.diff)[1]] <- 1
} 
H
}
