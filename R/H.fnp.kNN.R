H.fnp.kNN<-function(x,k=5,
semimetric="deriv",q=2,
kind.of.kernel="quad",range.grid=NULL,nknot=NULL)
{
if (!is.matrix(x)) x<-as.matrix(x)
p<-ncol(x)
kind.of.semimetric <- paste("semimetric.", semimetric, sep = "")
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
	else if (kind.of.semimetric=="semimetric.pca")  norm.diff <- semimetric.pca(data1=x[i,], data2=x, q=q)                      
	norm.order <- order(norm.diff[1,])
	ind.curves2 <- norm.order[1:(k + 1)]
	h <- sum(abs(norm.diff[ind.curves2[k:(k + 1)]]))/2
	res.kernel <- kernel(norm.diff[ind.curves2[ - (k + 1)]]/h)
	sum.res.kernel <- sum(res.kernel)
	for(j in 1:k){
		cur<-ind.curves2[j]
		if(sum.res.kernel > 0) {
			H[i,cur] <-res.kernel[j]/sum.res.kernel
		}                                                
		else H[i,cur] <- 1

	}
} 
H
}
