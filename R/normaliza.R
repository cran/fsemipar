normaliza<-function(coef,range.grid,t0=NULL,order.Bspline=3,nknot.theta=3)
{
if(is.vector(coef))coef <- as.matrix(t(coef))
if (ncol(coef) != (nknot.theta + order.Bspline)) stop("(length(theta) != (nknot.theta + order.Bspline))")
if (is.null(t0)) t0<-mean(range.grid)
a <- range.grid[1]
b <- range.grid[2]
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))     
point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142)
weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
x.gauss<- 0.5 * ((b + a) + (b - a) * point.gauss)
lx.gauss<- length(x.gauss)
Bspline.g.theta<-splineDesign(delta.theta, x.gauss, order.Bspline)
H <-t(Bspline.g.theta)%*%(Bspline.g.theta*(weight.gauss*0.5*(b-a)))
dim<-ncol(coef)
per<-nrow(coef)
theta.normal<-matrix(0,nrow=per,ncol=dim)
for(i in 1:per){
	coefi=coef[i,]
	prod.coef <-outer(coefi,coefi, "*")
	norm=sqrt(sum(H*prod.coef))
	theta.normal[i,]=coefi/norm
}
Bspline.t0.theta<-splineDesign(delta.theta, t0, order.Bspline)
theta.t0<-theta.normal%*%t(Bspline.t0.theta)
pos<-which(theta.t0>0) 
theta.normal.pos<-theta.normal[pos,]
theta.normal.pos
}


