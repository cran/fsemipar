H.fsim.kernel<-function(x,h,
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
	res.kernel <- kernel(norm.diff/h)
	res.kernel[res.kernel<0] <- 0
	res.kernel[res.kernel>1] <- 0
	sum.res.kernel <- sum(res.kernel)
	if(sum.res.kernel > 0) {
		H[i,] <-res.kernel/sum.res.kernel
	}                                                
	else H[i,which.min(norm.diff)] <- 1
} 
H
}



