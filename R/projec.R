projec<-function(data,
theta,order.Bspline=3, nknot.theta=3,
range.grid=NULL,nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if(is.vector(data)) data <- as.matrix(t(data))
p <- ncol(data)
if (is.null(range.grid)) range.grid <- c(1,p)
a <- range.grid[1]
b <- range.grid[2]
x <- seq(a, b, length = p)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
nknotmax <- (p - order.Bspline - 1)%/%2
if(nknot > nknotmax){
	stop(paste("Give a number nknot smaller than",nknotmax,"for avoiding ill-conditioned matrix."))
}
if(length(theta) != (nknot.theta + order.Bspline)) stop("(length(theta) != (nknot.theta + order.Bspline))")
Knot <-seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
delta <-sort(c(rep(c(a, b),order.Bspline), Knot))
Bspline <-splineDesign(delta,x,order.Bspline)
Cmat <-crossprod(Bspline)
Dmat1 <-crossprod(Bspline, t(data))
coef.mat1 <-symsolve(Cmat, Dmat1)
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))
Bspline.theta<-splineDesign(delta.theta,x,order.Bspline)
theta.rec<-crossprod(t(Bspline.theta),theta) 
Dmat.theta<-crossprod(Bspline,theta.rec)
Theta.coef<-symsolve(Cmat,Dmat.theta)
point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142)
weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
x.gauss<-0.5 * ((b + a) + (b - a) * point.gauss)
lx.gauss<-length(x.gauss)
Bspline.g<-splineDesign(delta, x.gauss, order.Bspline)
H <-t(Bspline.g)%*%(Bspline.g*(weight.gauss*0.5*(b-a)))
theta.x1<-matrix(0,nrow=length(Theta.coef),ncol=ncol(coef.mat1))
for(i in 1:ncol(coef.mat1)){
	theta.x1[,i]=coef.mat1[,i]*Theta.coef
}
coef <-t(H%*%theta.x1);coef
projec<-apply(coef,1,sum)
projec
}