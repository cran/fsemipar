fsim.kNN.fit<- function(x,y,
seed.coeff=c(-1,0,1), order.Bspline=3, nknot.theta=3, t0=NULL,
min.knn=2, max.knn=NULL, knearest=NULL,step=NULL,  
kind.of.kernel="quad",range.grid=NULL,nknot=NULL)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
if (is.null(max.knn)) max.knn <- (n %/% 2)
p <- ncol(x)
estimated.Y <- list()
length.curve.y<-ncol(y)
yhat.cv <- matrix(0,n,length.curve.y)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(t0)) t0 <- mean(range.grid)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
dim.base.theta <- order.Bspline + nknot.theta
THETA.seq <- permutations(length(seed.coeff), dim.base.theta, seed.coeff, repeats.allowed=TRUE)
THETA.seq <- THETA.seq[(apply(abs(THETA.seq) , 1,sum) != 0)  & (!is.na(apply(abs(THETA.seq), 1,sum))), ]
THETA.seq.normalizado <- normaliza(coef=THETA.seq, range.grid=range.grid, t0=t0, order.Bspline=order.Bspline, nknot.theta=nknot.theta)
num.norm <- nrow(THETA.seq.normalizado)
if (is.null(max.knn)) max.knn <- n%/%2
if (is.null(knearest)){
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
k.seq <- knearest
kmax <- max(k.seq)       
num.band <- length(k.seq)
k.vec <- rep(0,num.norm)                                                                        
CV.app <- 0
for(m in 1:num.norm){
	message(m,"/",num.norm)
	cv.kseq <- 0        
	for(i in 1:n) {                      
		norm.diff <- semimetric.projec(data1=x, data2=x[i,], theta=THETA.seq.normalizado[m,], range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)
		norm.order <- order(norm.diff)
		zz <- sort(norm.diff)[2:(kmax + 2)]
		bandwith <- 0.5 * (zz[-1] + zz[ - (kmax + 1)])
		z <- zz[ - (kmax + 1)]
		Zmat <- matrix(rep(z, kmax), nrow = kmax, byrow = T)
		Umat <- Zmat/bandwith
		Kmat <- kernel(Umat)
		Kmat[col(Kmat) > row(Kmat)] <- 0
		ind.curves1 <- norm.order[2:(kmax + 1)]
		yind <- y[ind.curves1,]
		if (length.curve.y==1) {                       
			Ymat <- matrix(rep(yind, kmax), nrow = kmax, byrow = T)
			yhat1 <- apply(Ymat[k.seq,  ] * Kmat[k.seq,  ], 1, sum)/apply(Kmat[k.seq,  ], 1, sum)
			resid.kseq.2 <- (yhat1- y[i])^2
			cv.kseq <- cv.kseq + resid.kseq.2 
		}
		else{
			Karr <- array(t(Kmat), c(dim(t(Kmat)),length.curve.y))
			Karr <- aperm(Karr, c(1,3,2))
			Yarr <- array(yind, c(dim(yind),kmax))
			yhat1 <- (t(apply(Yarr[,,k.seq] * Karr[,,k.seq], 2:3, sum)))/matrix(apply(Kmat[k.seq,], 1, sum),length(k.seq),length.curve.y)
			resid.kseq.2 <- apply((yhat1- matrix(rep(y[i,],nrow(yhat1)), dim(yhat1), byrow=TRUE))^2, 1, sum)  
			cv.kseq <- cv.kseq + resid.kseq.2
		}                  
	} 
	cv.kseq <- cv.kseq/(n*length.curve.y)
	index <- order(cv.kseq)[1]
	k.opt.m <- k.seq[index]
	k.vec[m] <- k.opt.m
	CV.app[m] <- cv.kseq[index]
	for (i in 1:n) 
		yhat.cv[i,] <- fsim.kNN.test(y=y[-i],x=x[-i,], y.test=y[i],x.test=x[i,], theta=THETA.seq.normalizado[m,], k=k.opt.m, kind.of.kernel=kind.of.kernel, range.grid=range.grid,
						order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)$y.estimated.test
		estimated.Y[[m]] <- yhat.cv
}
m.opt <- order(CV.app)[1]
theta.opt <- THETA.seq.normalizado[m.opt,]
k.opt <- k.vec[m.opt]
H<-H.fsim.kNN(x=x,k=k.opt,theta=theta.opt,nknot.theta=nknot.theta,order.Bspline=order.Bspline,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot)
CV.opt<-CV.app[m.opt]
fitted.values<-H%*%y
res<-y-drop(fitted.values)
df<-sum(diag(H))
sigma.res<-sum(res^2)/(n-df)
r.2<-1-sum(res^2)/sum((y-mean(y))^2)
call<-match.call()
out<-list(fitted.values=fitted.values,residuals=res,theta.est=theta.opt, k.opt=k.opt, r.squared=r.2,var.res=sigma.res, df=df,
		yhat.cv=estimated.Y[[order(CV.app)[1]]],CV.opt=CV.opt,CV.values=CV.app, H=H,
		m.opt=m.opt,theta.seq.norm=THETA.seq.normalizado, k.seq=k.seq,
		call=call,y=y,x=x,n=n,
		kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
		order.Bspline=order.Bspline,nknot.theta=nknot.theta)
class(out)<- "fsim.kNN"
out
}


