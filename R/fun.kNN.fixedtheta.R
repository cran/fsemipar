fun.kNN.fixedtheta <- function(y, x, pred, 
theta, order.Bspline=3, nknot.theta=3,
knearest=NULL, 
kind.of.kernel="quad",range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
y<- as.vector(y)
if (is.vector(x)) {
	x<- as.matrix(x)
	pred<- as.matrix(pred)
}
else if(is.vector(pred)) pred<- as.matrix(t(pred))
p <- ncol(x)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 		
kernel <- get(kind.of.kernel)
SEMIMETRIC1 <- semimetric.projec(data1=x, data2=pred, theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)
n1 <- ncol(SEMIMETRIC1)
step <- ceiling(n1/100)
if(step == 0) step <- 1
if (is.null(knearest)) knearest <- seq(from = 2, to = n1 %/% 2, by = step)
kmax <- max(knearest)		
SEMIMETRIC2 <- semimetric.projec(data1=x, data2=pred, theta=theta, range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)
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
	if (length(knearest)==1) Yhat[i,  ] <- sum(YMAT[knearest,  ] * KMAT[knearest,  ])/sum(KMAT[knearest,  ])
	else Yhat[i,  ] <- apply(YMAT[knearest,  ] * KMAT[knearest,  ],1, sum)/apply(KMAT[knearest,  ], 1, sum)
}	
list(yhat=Yhat, knn=knearest, h.seq=BANDWIDTH)
}















