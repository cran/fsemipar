fnp.kNN.fit.test.loc<- function(y,x,y.test=NULL,x.test,
kind.of.semimetric="semimetric.pca", start.order.deriv.o.pca=NULL, end.order.deriv.o.pca=NULL, min.leng.interv=NULL, max.leng.interv=NULL, 
knearest=NULL, min.knn=2, max.knn=NULL, step=NULL,
kind.of.kernel="quad", range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
kernel <- get(kind.of.kernel)
semimetric <- get(kind.of.semimetric)
if (!is.matrix(x))  x <- t(as.matrix(x))
if (!is.matrix(x.test))  x.test<- t(as.matrix(x.test))
n <- nrow(x)
p <- ncol(x)
J <- nrow(x.test)
if (is.null(max.knn)) max.knn <- n %/% 2
estimated.Y <- list()
if (p==1)  stop("A functional variable is needed (i.e., ncol(x)>2 is required")
if ((kind.of.semimetric=="semimetric.deriv") | (kind.of.semimetric=="semimetric.pca")){		       
	if ( (kind.of.semimetric=="semimetric.deriv") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 0
	if ( (kind.of.semimetric=="semimetric.deriv") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- 3
	if ((kind.of.semimetric=="semimetric.pca") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 1
	if ((kind.of.semimetric=="semimetric.pca") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- p
	start.order <- start.order.deriv.o.pca
	end.order <- end.order.deriv.o.pca
}		
else if (kind.of.semimetric=="semimetric.interv"){
	if (is.null(min.leng.interv)) min.leng.interv <- 1
	if (is.null(max.leng.interv)) max.leng.interv <- p-1
	start.order <- 1
	end.order <- (max.leng.interv-min.leng.interv +1)*(2*p-max.leng.interv-min.leng.interv)/2
}
else stop("kind.of.semimetric must be equal to: semimetric.deriv, semimetric.interv or semimetric.pca")
if (is.null(range.grid)) range.grid <- c(1,p)
nknot.m <- nknot
if (is.null(nknot)) {
	if (kind.of.semimetric=="semimetric.interv") nknot.m <- (p - 0 - 3 - 1)%/%2 
}
num.norm <- end.order - start.order +1
if (is.null(knearest)) {
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
k.seq <- knearest
kmax <- max(k.seq)
resid <- 0
CV.app <- 0
MSEP <- 0
k.mat <- matrix(0, n, num.norm)
for(m in start.order:end.order) {
	if (is.null(nknot)) {
		if (kind.of.semimetric=="semimetric.deriv") nknot.m <- (p - m - 3 - 1)%/%2 
	}
	a <- p-min.leng.interv
	if (kind.of.semimetric=="semimetric.interv") {
		for (i in 1:(max.leng.interv-min.leng.interv+1)) 
		if ( ((i-1)*(2*a+2-i)/2 + 1 <= m) & (m <= (i*(2*a+1-i)/2)) ) {		                  
			ext.inf <- m - (i-1)*(2*a+2-i)/2
			ext.sup <- ext.inf + min.leng.interv + i -1
		}
	}
	for(i in 1:n) {                      				
		if (kind.of.semimetric=="semimetric.deriv") norm.diff <- semimetric(data1=x, data2=x[i,], q=m, range.grid=range.grid, nknot = nknot.m)
		else if (kind.of.semimetric=="semimetric.interv") norm.diff <- semimetric(data1=x, data2=x[i,], interv=c(ext.inf,ext.sup), range.grid=range.grid, nknot = nknot.m)					  
		else if (kind.of.semimetric=="semimetric.pca")  norm.diff <- semimetric(data1=x, data2=x[i,], q=m)
		norm.order <- order(norm.diff)
		Y.obs1 <- y[i]
		zz <- sort(norm.diff)[2:(kmax + 2)]
		bandwith <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
		z <- zz[ - (kmax + 1)]
		Zmat <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
		Umat <- Zmat/bandwith
		Kmat <- kernel(Umat)
		Kmat[col(Kmat) > row(Kmat)] <- 0
		ind.curves1 <- norm.order[2:(kmax + 1)]
		yind <- y[ind.curves1]
		Ymat <- matrix(rep(yind, kmax), nrow = kmax, byrow = T)
		yhat1 <- rowSums(Ymat[k.seq,  ] * Kmat[k.seq,  ])/rowSums(Kmat[k.seq,  ])
		criterium <- abs(yhat1 - y[i])
		index <- order(criterium)[1]
		estimate <- yhat1[index]
		resid[i] <- estimate - Y.obs1
		k.mat[i, m + 1 - start.order] <- k.seq[index]
	} 			
	CV.app[m + 1 - start.order] <- sum(resid^2)/n
	y.hat2 <- 0
	for(j in 1:J) {												 
		if (kind.of.semimetric=="semimetric.deriv") norm.diff <- semimetric(data1=x, data2=x.test[j,], q=m, range.grid=range.grid, nknot = nknot.m)
		else if (kind.of.semimetric=="semimetric.interv") norm.diff <- semimetric(data1=x, data2=x.test[j,], interv=c(ext.inf,ext.sup), range.grid=range.grid, nknot = nknot.m)					  
		else if (kind.of.semimetric=="semimetric.pca")  norm.diff <- semimetric(data1=x, data2=x.test[j,], q=m)
		i.opt <- order(norm.diff)[1]
		k.opt <- k.mat[i.opt, m + 1 - start.order]
		norm.order <- order(norm.diff)
		ind.curves2 <- norm.order[1:(k.opt + 1)]
		bandwith.opt <- sum(abs(norm.diff[ind.curves2[k.opt:(k.opt + 1)]]))/2
		res.kernel <- kernel(norm.diff[ind.curves2[ - (k.opt + 1)]]/bandwith.opt)
		sum.res.kernel <- sum(res.kernel)
		if(sum.res.kernel > 0) y.hat2[j] <- sum(y[ind.curves2[ - (k.opt + 1)]] * res.kernel)/sum.res.kernel
		else y.hat2[j] <- y[ind.curves2[1]]
	}
	if(!is.null(y.test)){
		MSEP[m + 1 - start.order] <- sum((y.hat2 - y.test)^2)/J
	}
	else MSEP<-NULL
	estimated.Y[[m + 1 - start.order]] <- y.hat2	   
} # for (m
m.opt <- order(CV.app)[1] +  start.order - 1
if (kind.of.semimetric=="semimetric.interv") {		
	a <- p-min.leng.interv   	
	for (i in 1:(max.leng.interv-min.leng.interv+1)) 
		if ( ((i-1)*(2*a+2-i)/2 + 1 <= m.opt) & (m.opt <= (i*(2*a+1-i)/2)) ) {
			ext.inf <- m.opt - (i-1)*(2*a+2-i)/2
			ext.sup <- ext.inf + min.leng.interv + i -1
		}
	list(m.opt=m.opt, interv.opt=c(ext.inf, ext.sup), MSE.test.m.opt=MSEP[order(CV.app)[1]], y.estimate.test = estimated.Y[[order(CV.app)[1]]], 
	CV.app=CV.app, MSE.test = MSEP, y.test=y.test)
}
else
	list(m.opt=m.opt, MSE.test.m.opt=MSEP[order(CV.app)[1]], y.estimate.test = estimated.Y[[order(CV.app)[1]]], 
	CV.app=CV.app, MSE.test = MSEP,y.test=y.test)	
}





