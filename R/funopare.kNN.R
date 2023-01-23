funopare.kNN <- function(y, x, pred, k, ... , kind.of.kernel = "quad", semimetric = "pca")
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
y <- as.vector(y)
if(is.vector(pred)) pred <- as.matrix(t(pred))
testfordim <- sum(dim(x)==dim(pred))==2
twodatasets <- T
if(testfordim) twodatasets <- sum(x==pred)!=prod(dim(x))
sm <- get(paste("semimetric.", semimetric, sep = ""))
if(semimetric == "mplsr")
	SEMIMETRIC1 <- sm(y, x, x, ...)
else SEMIMETRIC1 <- sm(x, x, ...)
kernel <- get(kind.of.kernel)
p1 <- ncol(SEMIMETRIC1)
n1 <- nrow(SEMIMETRIC1)
if(k >= n1)
	stop(paste("try a smaller number of neighbour \n than ", k))
bandwidth.knn1 <- 0
for(j in 1:p1) {
	Sem <- SEMIMETRIC1[, j]
	knn.to.band <- Sem[order(Sem)[k:(k + 1)]]
	bandwidth.knn1[j] <- 0.5 * sum(knn.to.band)
}
KERNEL1 <- kernel(t(t(SEMIMETRIC1)/bandwidth.knn1))
KERNEL1[KERNEL1 < 0] <- 0
KERNEL1[KERNEL1 > 1] <- 0
diag(KERNEL1) <- 0
RESPKERNEL1 <- KERNEL1 * y
Denom1 <- apply(KERNEL1, 2, sum)
Response.estimated <- apply(RESPKERNEL1, 2, sum)/Denom1
Mse.estimated <- sum((Response.estimated - y)^2)/n1
if(twodatasets) {
	if(semimetric == "mplsr")
		SEMIMETRIC2 <- sm(y, x, pred, ...)
	else SEMIMETRIC2 <- sm(x, pred, ...)
	p2 <- ncol(SEMIMETRIC2)
	bandwidth.knn2 <- 0
	for(j in 1:p2) {
		Sem <- SEMIMETRIC2[, j]
		knn.to.band <- Sem[order(Sem)[k:(k + 1)
			]]
		bandwidth.knn2[j] <- 0.5 * sum(knn.to.band)
	}
	KERNEL2 <- kernel(t(t(SEMIMETRIC2)/bandwidth.knn2))
	KERNEL2[KERNEL2 < 0] <- 0
	KERNEL2[KERNEL2 > 1] <- 0
	Denom2 <- apply(KERNEL2, 2, sum)
	RESPKERNEL2 <- KERNEL2 * y
	Response.predicted <- apply(RESPKERNEL2, 2, sum)/Denom2
	return(list(estimated.values = Response.estimated, 
		predicted.values = Response.predicted, knn = k, 
		mse = Mse.estimated))
}
else {
	return(list(estimated.values = Response.estimated, knn = 
		k, mse = Mse.estimated))
}
}
