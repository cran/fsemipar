H.fsim.kNN<-function(x,k,
theta,nknot.theta=3,order.Bspline=3,
kind.of.kernel="quad",range.grid=NULL,nknot=NULL)
{
if (!is.matrix(x)) x<-as.matrix(x)
p<-ncol(x)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
n<-nrow(x)
kernel <- get(kind.of.kernel)
H<-matrix(0,nrow=n,ncol=n)
for(i in 1:n){
	norm.diff <- semimetric.projec(data1=x[i,], data2=x, theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot = nknot, nknot.theta=nknot.theta)
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




