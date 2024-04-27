fun.kNN.fixedtheta<- function(y, norm.diff,knearest=NULL,kernel=kernel)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
y<- as.vector(y)	
n2<-ncol(norm.diff)
step <- ceiling(n2/100)
if(step == 0) step <- 1
if (is.null(knearest)) knearest <- seq(from = 2, to = n2 %/% 2, by = step)
kmax <- max(knearest)		
Yhat <- matrix(0, nrow = n2, ncol = length(knearest))
BANDWIDTH <- matrix(0, nrow = n2, ncol = kmax)
 
for(i in 1:n2) {
	Norm.diff <- norm.diff[, i]		
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
	Yhat[i,  ] <- rowSums(YMAT[knearest,  ] * KMAT[knearest,  ])/rowSums(KMAT[knearest,  ])
}	
list(yhat=Yhat, knn=knearest, h.seq=BANDWIDTH)
}















