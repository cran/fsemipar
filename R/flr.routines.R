#############################################################################
#############################################################################
#############################################################################
### FUNCTIONAL LINEAR REGRESSION (see Cardot, Ferraty and Sarda, 2003)   ####
### The routines given just below translate the S+ routines downloadable ####
### at http://www.math.univ-toulouse.fr/staph/ENGLISH/a-logiciels.html   ####
### (file "pgsplus.spl") into R language.                                ####
#############################################################################
#############################################################################
#############################################################################
Splinemlf <- function(Y, X, lambda, Posintknots, order = 4, m = 2)
{
n <- nrow(X)
p <- ncol(X)
res <- Bspline.ini(t(X), Posintknots, order, m)
Gamma <- res$A %*% t(res$A)/n + lambda * res$G
Delta <- t(as.matrix(Y)) %*% t(res$A)/n
alpha0 <- Delta %*% solve(Gamma)
alpha1 <- res$B %*% t(alpha0)
yhat <- as.vector(t(alpha1) %*% t(X)/p)
list(as = alpha0, af = alpha1, yhat = yhat, rd = Y - yhat, r2 = 1 - mean((Y - yhat)^2)/var(Y))
}


Bspline.ini <- function(X, Posintknots, order, m)
{
p <- nrow(X)
n <- ncol(X)
nknot <- length(Posintknots)
A <- matrix(ncol = n, nrow = order + nknot)
x0 <- seq(0, 1, length = p)
x <- seq(0, 1, length = 200)
Intknots <- x0[Posintknots]
delta <- sort(c(rep(range(x), order), Intknots))
B <- splineDesign(delta, x, order)
xdiff <- diff(x, 1)
DmBj <- splineDesign(delta, x, order, derivs = rep(m, 200))
G1 <- t(DmBj[-1,  ]) %*% (DmBj[-1,  ] * xdiff)
G2 <- t(DmBj[-200,  ]) %*% (DmBj[-200,  ] * xdiff)
G <- 0.5 * (G1 + G2)    ## Calcul des coordonnees
B <- splineDesign(delta, x0, order)
A <- t(B) %*% X/p
list(A = A, G = G, B = B)
}


approx.spline.deriv <- function(DATA, nderivs, nknot, Range.grid, degree=3)
{
p <- ncol(DATA)
a <- Range.grid[1]
b <- Range.grid[2]
Knot <- seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]
x <- seq(a, b, length = p)
ord <- nderivs + degree
nknotmax <- (p - ord - 1)%/%2
if(nknot > nknotmax){
	stop(paste("give a number nknot smaller than ",nknotmax, " for avoiding ill-conditioned matrix"))
}
Delta <- sort(c(rep(c(a, b), ord), Knot))
BSPLINE <- splineDesign(Delta, x, ord)
CMAT <- crossprod(BSPLINE)
DMAT <- crossprod(BSPLINE, t(DATA))
COEF <- symsolve(CMAT, DMAT)
BSPLINEDER <- splineDesign(Delta, x, ord, rep(nderivs, length(x)))
list(COEF = COEF, APPROX = BSPLINEDER %*% COEF)
}


interp.spline.deriv <- function(DATA, nderivs, intknot, Design1, Design2, degree=3)
{
p <- ncol(DATA)
a <- min(Design1)
b <- max(Design1)
if(length(intknot)==1){
	Knot <- seq(a, b, length = intknot + 2)[ - c(1, intknot + 2)]
}
else{     
	Knot <- intknot
}
ord <- nderivs + degree
nknotmax <- (p - ord - 1)%/%2
Delta <- sort(c(rep(c(a, b), ord), Knot))
BSPLINE <- splineDesign(Delta, Design2, ord)
CMAT <- crossprod(BSPLINE)
DMAT <- crossprod(BSPLINE, t(DATA))
COEF <- symsolve(CMAT, DMAT)
BSPLINEDER <- splineDesign(Delta, Design1, ord, rep(nderivs, length(Design1)))
list(COEF = COEF, APPROX = BSPLINEDER %*% COEF)
}
