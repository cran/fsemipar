IASSMR.kernel.fit<- function(x, z, y, train.1=NULL, train.2=NULL, 
seed.coeff=c(-1,0,1), order.Bspline=3, nknot.theta=3, t0=NULL, 
min.q.h=0.05, max.q.h=0.5, h.seq = NULL, num.h = 10,
range.grid=NULL, kind.of.kernel="quad", nknot=NULL,
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, vn=ncol(z), nfolds=10, seed=123, wn=c(10,15,20), 
criterion=c("GCV", "BIC", "AIC", "k-fold-CV"),
penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"), 
max.iter=1000)
{
if (!is.matrix(z)) z <- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
n<-nrow(z)
min.quantile.h<-min.q.h
max.quantile.h<-max.q.h
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
if(is.null(train.1)) train.1<-1:ceiling(n/2)
if(is.null(train.2)) train.2<-(ceiling(n/2)+1):n
n1 <- length(train.1)
n2 <- length(train.2)
num.wn <- length(wn)
pn <- ncol(z)
indexes.beta <- 1:pn
p<-ncol(x)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
if (!(is.null(h.seq))) num.h <- length(h.seq)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(t0)) t0 <- mean(range.grid)
index01 <- list()
index1 <- list()
lambda1 <- 0
h1 <- 0
beta1 <- list()
theta1 <- list()
IC1 <- rep(Inf, length=num.wn)
vn1<-numeric(num.wn)
index2 <- list()
lambda2 <- 0
h2 <- 0
beta2 <- list()
theta2 <- list()
IC2 <- rep(Inf, length=num.wn)
vn2<-numeric(num.wn)
indexes.beta.nonnull2 <- list()
prueba<-list()
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
	step.1 <-sfplsim.kernel.fit(x=x[train.1,],z=z[train.1,index.1], y=y[train.1], 
								seed.coeff=seed.coeff, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta, t0=t0, 
								lambda.min=lambda.min, lambda.min.h=lambda.min.pn.high, lambda.min.l=lambda.min.pn.low, factor.pn=factor.pn,
								nlambda=nlambda, lambda.seq=NULL, vn=vn, nfolds=nfolds, seed=seed, 
								min.q.h=min.quantile.h, max.q.h=max.quantile.h, h.seq = h.seq, num.h = num.h,
								range.grid=range.grid, kind.of.kernel=kind.of.kernel,
								criterion=criterion, penalty=penalty, 
								max.iter=max.iter)
	beta <- step.1$beta.est
	index.X.pen <- step.1$indexes.beta.nonnull
	beta1[[w]] <- beta
	theta1[[w]] <- step.1$theta.est
	index1[[w]] <- index.1[index.X.pen]
	lambda1[w] <- step.1$lambda.opt 
      vn1[w]<-step.1$vn.opt
	h1[w] <- step.1$h.opt
	IC1[w] <- step.1$IC
	if  (sum(index.X.pen)==0){		
		beta2[[w]] <- NaN
		indexes.beta.nonnull2[[w]] <- NaN
		IC2[w] <- NaN
		lambda2[w] <- NaN
		h2[w] <- NaN
		index2[[w]] <- NaN
		next
	}		
	index.2 <- 0
	index.infl <- (1:wn[w])[index.X.pen]
	for (j in index.infl) index.2 <- c(index.2, index.wn[[j]])
	index.2 <- index.2[-1]
	index2[[w]] <- index.2
	step.2 <- sfplsim.kernel.fit(x=x[train.2,],z=z[train.2,index.2], y=y[train.2],
								 seed.coeff=seed.coeff, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta, t0=t0, 
								 lambda.min=lambda.min, lambda.min.h=lambda.min.pn.high, lambda.min.l=lambda.min.pn.low, factor.pn=factor.pn,
								 nlambda=nlambda, lambda.seq=NULL, vn=vn, nfolds=nfolds, seed=seed, 
								 min.q.h=min.quantile.h, max.q.h=max.quantile.h, h.seq = h.seq, num.h = num.h,
								 range.grid=range.grid, kind.of.kernel=kind.of.kernel,
								 criterion=criterion, penalty=penalty, 
								 max.iter=max.iter)
	beta2[[w]] <- step.2$beta.est
	theta2[[w]] <- step.2$theta.est
	indexes.beta.nonnull2[[w]] <- indexes.beta[index.2][step.2$indexes.beta.nonnull]
	lambda2[w] <- step.2$lambda.opt
	h2[w] <- step.2$h.opt
	IC2[w] <- step.2$IC
      vn2[w]<-step.2$vn.opt			
} 
index.w.opt<- order(IC2)[1]
beta.red<-beta2[[index.w.opt]]
w.opt<-wn[index.w.opt]
theta.est<-theta2[[index.w.opt]]
IC<-IC2[index.w.opt]
vn.opt<-vn2[index.w.opt]
lambda.opt<-lambda2[index.w.opt]
h.opt<-h2[index.w.opt]
indexes.beta.nonnull<-indexes.beta.nonnull2[[index.w.opt]]
beta.red.nonnull<-beta.red[beta.red!=0]
beta.est<-numeric(pn)
beta.est[indexes.beta.nonnull]<-beta.red.nonnull
ww<-H.fsim.kernel(x=x,h=h.opt,theta=theta.est,nknot.theta=nknot.theta,order.Bspline=order.Bspline,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot)
yhp<-z%*%beta.est+ww%*%(y-z%*%beta.est)
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta.est, theta.est=theta.est,indexes.beta.nonnull=indexes.beta.nonnull,w.opt=w.opt,h.opt=h.opt,
		  lambda.opt=lambda.opt,IC=IC,vn.opt=vn.opt,
		  beta2=beta2,theta2=theta2, indexes.beta.nonnull2=indexes.beta.nonnull2, h2=h2, IC2=IC2,lambda2=lambda2,index02=index2, 
		  beta1=beta1,theta1=theta1,h1=h1, IC1=IC1, lambda1=lambda1, index01=index01, index1=index1,
		  call=call,y=y,x=x,z=z,train.1=train.1,train.2=train.2,
		  kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
		  order.Bspline=order.Bspline,nknot.theta=nknot.theta,seed.coeff=seed.coeff,t0=t0,
		  criterion=criterion, penalty=penalty,max.iter=max.iter,wn=wn,
		  seed=seed,nlambda=nlambda,lambda.min.pn.low=lambda.min.pn.low,lambda.min.pn.high=lambda.min.pn.high,factor.pn=factor.pn, 
		  nfolds=nfolds,group=group,vn=vn,
		  num.h=num.h,h.seq=h.seq,min.quantile.h=min.quantile.h,max.quantile.h=max.quantile.h)
class(out)<-"IASSMR.kernel"
out
}

