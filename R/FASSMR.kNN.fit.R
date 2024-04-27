FASSMR.kNN.fit<- function(x,z,y,  
seed.coeff=c(-1,0,1), order.Bspline=3,  nknot.theta=3, 
knearest=NULL, min.knn=2, max.knn=NULL, step=NULL,
kind.of.kernel="quad",range.grid=NULL, nknot=NULL, 
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, vn=ncol(z), nfolds=10, seed=123, wn=c(10,15,20),
criterion="GCV",
penalty="grSCAD", 
max.iter=1000,n.core=NULL)
{
if (!is.matrix(z)) z<- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
n<-nrow(x)
num.wn <- length(wn)
pn <- ncol(z)
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
p<-ncol(x)
indexes.beta <- 1:pn
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
if (is.null(range.grid)) range.grid <- c(1,p)
index01 <- vector("list",length=num.wn)
index1 <- vector("list",length=num.wn)
lambda1 <- 0
knn1 <- 0
beta1 <- vector("list",length=num.wn)
theta1 <- vector("list",length=num.wn)
ww<-vector("list",length=num.wn)
norm.diff<-vector("list",length=num.wn)
IC1 <- rep(Inf, length=num.wn)
indexes.beta.nonnull1 <- vector("list",length=num.wn)
vn1<-numeric(num.wn)
for (w in 1:num.wn) {
	message("wn=", wn[w], ": ", w, "/", num.wn)
	num.veci <- trunc(pn/wn[w])
	aux <- pn - wn[w]*num.veci
	group <- 0
	if  (aux!=0) {
		for (j in 1:aux) group <- c(group, rep(j, length=num.veci+1))
		for (j in (aux+1):wn[w]) group <- c(group, rep(j, length=num.veci))
	}
	else for (j in 1:wn[w]) group <- c(group, rep(j, length=num.veci))
	group <- group[-1]
	index.wn <- list()
	index.1 <- 0
	for (j in 1:wn[w]) {	
		index.wn[[j]] <- (1:pn)[group==j]
		index.1 <- c(index.1, trunc(median(index.wn[[j]])))
	}
	index.1 <- index.1[-1]
	index01[[w]] <- index.1
	step.1 <- sfplsim.kNN.fit(x=x, z=z[,index.1], y=y,  
							  seed.coeff=seed.coeff, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta, 
							  lambda.min=lambda.min, lambda.min.h=lambda.min.pn.high, lambda.min.l=lambda.min.pn.low, factor.pn=factor.pn,
							  nlambda=nlambda, lambda.seq=NULL, vn=vn, nfolds=nfolds, seed=seed, 
							  knearest=knearest, min.knn=min.knn, max.knn=max.knn, step=step,
							  range.grid=range.grid, kind.of.kernel=kind.of.kernel,
							  criterion=criterion, penalty=penalty, 
							  max.iter=max.iter,n.core=n.core)
	beta <- step.1$beta.est
	index.X.pen <- step.1$indexes.beta.nonnull
	beta1[[w]] <- beta
	theta1[[w]] <- step.1$theta.est
	index1[[w]] <- index.1[index.X.pen]
	lambda1[w] <- step.1$lambda.opt 
	knn1[w] <- step.1$k.opt
	IC1[w] <- step.1$IC
    vn1[w]<-step.1$vn.opt
	indexes.beta.nonnull1[[w]]<-index1[[w]]
	ww[[w]]<-step.1$ww
	norm.diff[[w]]<-step.1$norm.diff	
} 
index.w.opt<- which.min(IC1)
beta.red<-beta1[[index.w.opt]]
vn.opt<-vn1[index.w.opt]
w.opt<-wn[index.w.opt]
theta.est<-theta1[[index.w.opt]]
IC<-IC1[index.w.opt]
lambda.opt<-lambda1[index.w.opt]
k.opt<-knn1[index.w.opt]
indexes.beta.nonnull<-indexes.beta.nonnull1[[index.w.opt]]
beta.red.nonnull<-beta.red[beta.red!=0]
beta.est<-numeric(pn)
beta.est[indexes.beta.nonnull]<-beta.red.nonnull
ww<-ww[[index.w.opt]]
norm.diff<-norm.diff[[index.w.opt]]
yhp<-z%*%beta.est+ww%*%(y-z%*%beta.est)
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta.est,beta.red=beta.red,theta.est=theta.est,indexes.beta.nonnull=indexes.beta.nonnull,k.opt=k.opt,w.opt=w.opt,lambda.opt=lambda.opt,IC=IC,vn.opt=vn.opt,
		  beta.w=beta1,theta.w=theta1, IC.w=IC1,indexes.beta.nonnull.w=indexes.beta.nonnull1,lambda.w=lambda1, k.w=knn1, index01=index01,norm.diff=norm.diff, 
		  call=call,y=y,x=x,z=z,
		  kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
		  order.Bspline=order.Bspline,nknot.theta=nknot.theta,seed.coeff=seed.coeff,
		  criterion=criterion, penalty=penalty,max.iter=max.iter, wn=wn,
		  seed=seed,nlambda=nlambda,lambda.min.pn.low=lambda.min.pn.low, lambda.min.pn.high=lambda.min.pn.high, factor.pn=factor.pn,
		  nfolds=nfolds,vn=vn,
		  step=step,knearest=knearest,min.knn=min.knn,max.knn=max.knn)
class(out)<-"FASSMR.kNN"
out
}