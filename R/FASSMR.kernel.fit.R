FASSMR.kernel.fit<- function(x,z,y, 
seed.coeff=c(-1,0,1), order.Bspline=3,  nknot.theta=3, t0=NULL,
min.q.h=0.05, max.q.h=0.5, h.seq = NULL, num.h = 10,
range.grid=NULL, kind.of.kernel="quad", nknot=NULL,
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, vn=ncol(z), nfolds=10, seed=123,  wn=c(10,15,20),
criterion=c("GCV", "BIC", "AIC", "k-fold-CV"),  
penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso",
"gMCP"), max.iter=1000)
{
if (!is.matrix(z)) z<- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
n<-nrow(x)
num.wn <- length(wn)
pn <- ncol(z)
min.quantile.h<-min.q.h
max.quantile.h<-max.q.h
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
p<-ncol(x)
indexes.beta <- 1:pn
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
indexes.beta.nonnull1 <- list()
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
	step.1 <-sfplsim.kernel.fit(x=x, z=z[,index.1], y=y,  
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
	h1[w] <- step.1$h.opt
	IC1[w] <- step.1$IC
      vn1[w]<-step.1$vn.opt
	indexes.beta.nonnull1[[w]]<-index1[[w]]				
} 
index.w.opt<- order(IC1)[1]
beta.red<-beta1[[index.w.opt]]
w.opt<-wn[index.w.opt]
vn.opt<-vn1[index.w.opt]
theta.est<-theta1[[index.w.opt]]
IC<-IC1[index.w.opt]
lambda.opt<-lambda1[index.w.opt]
h.opt<-h1[index.w.opt]
indexes.beta.nonnull<-indexes.beta.nonnull1[[index.w.opt]]
beta.red.nonnull<-beta.red[beta.red!=0]
beta.est<-numeric(pn)
beta.est[indexes.beta.nonnull]<-beta.red.nonnull
ww<-H.fsim.kernel(x=x,h=h.opt,theta=theta.est,nknot.theta=nknot.theta,order.Bspline=order.Bspline,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot)
yhp<-z%*%beta.est+ww%*%(y-z%*%beta.est)
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta.est,beta.red=beta.red, theta.est=theta.est, indexes.beta.nonnull=indexes.beta.nonnull,h.opt=h.opt,
		w.opt=w.opt,lambda.opt=lambda.opt,IC=IC,vn.opt=vn.opt,
		beta.w=beta1,theta.w=theta1, IC.w=IC1,indexes.beta.nonnull.w=indexes.beta.nonnull1,lambda.w=lambda1,h.w=h1,index01=index01,
		call=call,y=y,x=x,z=z,
		kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
		order.Bspline=order.Bspline,nknot=nknot,nknot.theta=nknot.theta,seed.coeff=seed.coeff,t0=t0,
		criterion=criterion, penalty=penalty, max.iter=max.iter,wn=wn,
		seed=seed,nlambda=nlambda,lambda.min.pn.low=lambda.min.pn.low,lambda.min.pn.high=lambda.min.pn.high,factor.pn=factor.pn, 
		nfolds=nfolds,vn=vn,
		num.h=num.h,h.seq=h.seq,min.quantile.h=min.quantile.h,max.quantile.h=max.quantile.h)
class(out)<-"FASSMR.kernel"
out
}
