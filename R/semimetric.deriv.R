semimetric.deriv <- function(data1, data2, 
q=0, range.grid=c(1,20), nknot=10)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if(is.vector(data1)) data1 <- as.matrix(t(data1))
if(is.vector(data2)) data2 <- as.matrix(t(data2))
testfordim <- sum(dim(data1)==dim(data2))==2
twodatasets <- T
if(testfordim) twodatasets <- sum(data1==data2)!=prod(dim(data1))
p <- ncol(data1)
a <- range.grid[1]
b <- range.grid[2]
x <- seq(a, b, length = p)
order.Bspline <- q + 3
nknotmax <- (p - order.Bspline - 1)%/%2
if(nknot > nknotmax){
	stop(paste("give a number nknot smaller than ",nknotmax, " for avoiding ill-conditioned matrix"))
}
Knot <- seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
delta <- sort(c(rep(c(a, b), order.Bspline), Knot))
Bspline <- splineDesign(delta, x, order.Bspline)
Cmat <- crossprod(Bspline)
Dmat1 <- crossprod(Bspline, t(data1))
coef.mat1 <- symsolve(Cmat, Dmat1)
point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142)
weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
x.gauss <- 0.5 * ((b + a) + (b - a) * point.gauss)
lx.gauss <- length(x.gauss)
Bspline.deriv <- splineDesign(delta, x.gauss, order.Bspline, rep(q, lx.gauss))
H <- t(Bspline.deriv) %*% (Bspline.deriv * (weight.gauss * 0.5 * (b - a)))
eigH <- eigen(H, symmetric = T)
eigH$values[eigH$values < 0] <- 0
Hhalf <- t(eigH$vectors %*% (t(eigH$vectors) * sqrt(eigH$values)))
COEF1 <- t(Hhalf %*% coef.mat1)
if(twodatasets){
	Dmat2 <- crossprod(Bspline, t(data2))
	coef.mat2 <- symsolve(Cmat, Dmat2)
	COEF2 <- t(Hhalf %*% coef.mat2)
}
else{
	COEF2 <- COEF1
}
SEMIMETRIC <- 0
nbasis <- nrow(H)
for(f in 1:nbasis) SEMIMETRIC <- SEMIMETRIC + outer(COEF1[, f], COEF2[, f], "-")^2
out<-sqrt(SEMIMETRIC)
out
}
