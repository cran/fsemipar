semimetric.pca <- function(data1, data2, q=10)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if(is.vector(data1)) data1 <- as.matrix(t(data1))
if(is.vector(data2)) data2 <- as.matrix(t(data2))
testfordim <- sum(dim(data1)==dim(data2))==2
twodatasets <- T
if(testfordim) twodatasets <- sum(data1==data2)!=prod(dim(data1))
qmax <- ncol(data1)
if(q > qmax) stop(paste("give a integer q smaller than ", qmax))
n <- nrow(data1)
COVARIANCE <- t(data1) %*% data1/n
EIGENVECTORS <- eigen(COVARIANCE, symmetric = T)$vectors[,1:q]
COMPONENT1 <- data1 %*% EIGENVECTORS
if(twodatasets) {
	COMPONENT2 <- data2 %*% EIGENVECTORS
}
else {
	COMPONENT2 <- COMPONENT1
}
SEMIMETRIC <- 0
for(qq in 1:q)
	SEMIMETRIC <- SEMIMETRIC + outer(COMPONENT1[, qq], COMPONENT2[, qq], "-")^2
return(sqrt(SEMIMETRIC))
}