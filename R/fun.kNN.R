fun.kNN <- function(y, x, pred, 
semimetric = "deriv", q=2, 
knearest=NULL, 
range.grid=NULL,  kind.of.kernel="quad",nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
y<- as.vector(y)
if (is.vector(x)){
	x <- as.matrix(x)
	pred <- as.matrix(pred)
}
else if(is.vector(pred)) pred <- as.matrix(t(pred))
p <- ncol(x)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot)) nknot <- (p - q - 3 - 1)%/%2 
kernel <- get(kind.of.kernel)
if (p>1){
	sm <- get(paste("semimetric.", semimetric, sep = ""))
	if (semimetric=="deriv") 
	SEMIMETRIC1 <- sm(x, x, q=q, range.grid=range.grid, nknot=nknot)
	else SEMIMETRIC1 <- sm(x, x, q=q)
}
else{
	SEMIMETRIC1 <- matrix(0,nrow(x),nrow(x))
	for (i in 1:nrow(x)) SEMIMETRIC1[,i] <- abs(x-x[i])
}
n1 <- ncol(SEMIMETRIC1)
step <- ceiling(n1/100)
if(step == 0) step <- 1
if (is.null(knearest)) knearest <- seq(from = 2, to = n1 %/% 2, by = step)
kmax <- max(knearest)	
if (p>1) {
	if (semimetric=="deriv") SEMIMETRIC2 <- sm(x, pred, q=q, range.grid=range.grid, nknot=nknot)
	else SEMIMETRIC2 <- sm(x, pred, q=q)
}
else{
	SEMIMETRIC2 <- matrix(0,nrow(x),nrow(pred))
	for (i in 1:nrow(pred)) SEMIMETRIC2[,i] <- abs(x-pred[i])
}
n2 <- ncol(SEMIMETRIC2)
Yhat <- matrix(0, nrow = n2, ncol = length(knearest))
BANDWIDTH <- matrix(0, nrow = n2, ncol = kmax)
for(i in 1:n2) {
	Norm.diff <- SEMIMETRIC2[, i]	
	Norm.order <- order(Norm.diff)	
	zz <- sort(Norm.diff)[2:(kmax + 2)]	
	BANDWIDTH[i,  ] <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
	z <- zz[ - (kmax + 1)]
	ZMAT <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
	UMAT <- ZMAT/BANDWIDTH[i,  ]
	KMAT <- kernel(UMAT)
	KMAT[col(KMAT) > row(KMAT)] <- 0
	Ind.curves <- Norm.order[2:(kmax + 1)]
	Ind.resp <- y[Ind.curves]
	YMAT <- matrix(rep(Ind.resp, kmax), nrow = kmax, byrow = T)
	if (length(knearest)==1) Yhat[i,] <- sum(YMAT[knearest,] * KMAT[knearest,])/sum(KMAT[knearest,])
	else Yhat[i,] <- rowSums(YMAT[knearest,] * KMAT[knearest,])/rowSums(KMAT[knearest,])
}	
list(yhat=Yhat, knn=knearest, h.seq=BANDWIDTH)
}


















