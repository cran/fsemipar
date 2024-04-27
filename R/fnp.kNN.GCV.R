fnp.kNN.GCV<- function(y, x, pred, 
semimetric = "deriv",
knearest=NULL, min.knn=NULL, max.knn=NULL, step=NULL, ..., 
kind.of.kernel = "quad")
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
y <- as.vector(y)
if(is.vector(pred)) pred <- as.matrix(t(pred))
testfordim <- sum(dim(x)==dim(pred))==2
twodatasets <- T
if(testfordim) twodatasets <- sum(x==pred)!=prod(dim(x))
sm <- get(paste("semimetric.", semimetric, sep = ""))
if(semimetric == "mplsr") SEMIMETRIC1 <- sm(y, x, x, ...)
else SEMIMETRIC1 <- sm(x, x, ...)
kernel <- get(kind.of.kernel)
n1 <- ncol(SEMIMETRIC1)
if (is.null(min.knn)) min.knn <- 2
if (is.null(max.knn)) max.knn <- n1%/%2
if (is.null(knearest)) {
	if (is.null(step)) step <- ceiling(n1/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
kmax <- max(knearest)	
y.estimated <- 0
Bandwidth.opt <- 0
HAT.RESP <- matrix(0, nrow = n1, ncol = length(knearest))
BANDWIDTH <- matrix(0, nrow = n1, ncol = kmax)
for(i in 1:n1) {
	Norm.diff <- SEMIMETRIC1[, i]	
	Norm.order <- order(Norm.diff)	
	zz <- sort(Norm.diff)[2:(kmax + 2)]	
	BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
	z <- zz[ - (kmax + 1)]
	ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
	UMAT <- ZMAT/BANDWIDTH[i,  ]
	KMAT <- kernel(UMAT)
	KMAT[col(KMAT) > row(KMAT)] <- 0
	Ind.x <- Norm.order[2:(kmax + 1)]
	Ind.resp <- y[Ind.x]
	YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
	HAT.RESP[i,] <- rowSums(YMAT[knearest,] * KMAT[knearest,])/rowSums(KMAT[knearest,])
}
CRITERIUM <- (HAT.RESP - y)^2
Criterium <- colSums(CRITERIUM)
index.opt <- order(Criterium)[1]
y.estimated <- HAT.RESP[, index.opt]
knearest.opt <- knearest[index.opt]
Bandwidth.opt <- BANDWIDTH[, knearest.opt]
Mse.estimated <- sum((y.estimated - y)^2)/n1
if(twodatasets) {
	if(semimetric == "mplsr") SEMIMETRIC2 <- sm(y, x, pred, ...)
	else SEMIMETRIC2 <- sm(x, pred, ...)
	Bandwidth2 <- 0
	n2 <- ncol(SEMIMETRIC2)
	for(k in 1:n2) {
		Sm2k <- SEMIMETRIC2[, k]
		Bandwidth2[k] <- sum(sort(Sm2k)[knearest.opt:(knearest.opt+1)])*0.5
	}
	KERNEL <- kernel(t(t(SEMIMETRIC2)/Bandwidth2))
	KERNEL[KERNEL < 0] <- 0
	KERNEL[KERNEL > 1] <- 0
	Denom <- colSums(KERNEL)
	RESPKERNEL <- KERNEL * y
	y.predicted <- colSums(RESPKERNEL)/Denom
	return(list(Estimated.values = y.estimated, 	
		predicted.values = y.predicted, Bandwidths = 
		Bandwidth.opt, knearest.opt = knearest.opt, Mse = 
		Mse.estimated))
}
else {
	return(list(Estimated.values = y.estimated, Bandwidths
		= Bandwidth.opt, knearest.opt = knearest.opt, Mse = 
		Mse.estimated))
}
}
