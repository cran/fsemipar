PVS.kNN.fit<- function(x, z, y, train.1=NULL, train.2=NULL, 
semimetric="deriv", q=NULL,
knearest=NULL, min.knn=2, max.knn=NULL, step=NULL,
range.grid=NULL, kind.of.kernel="quad",nknot=NULL,
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, vn=ncol(z), nfolds=10, seed=123, wn=c(10,15,20),
criterion=c("GCV", "BIC", "AIC", "k-fold-CV"), 
penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"), 
max.iter=1000)
{
if (is.null(semimetric)) semimetric <- "deriv"
if (!is.matrix(z)) z <- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
kind.of.semimetric <- paste("semimetric.", semimetric, sep = "")
if (is.null(q)) q<-ifelse(semimetric=="deriv",0,2)
n<-nrow(z)
if(is.null(train.1)) train.1<-1:ceiling(n/2)
if(is.null(train.2)) train.2<-(ceiling(n/2)+1):n
n1 <- length(train.1)
n2 <- length(train.2)
num.wn <- length(wn)
pn <- ncol(z)
indexes.beta <- 1:pn
p<-ncol(x)
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
knn.opt <- 0
index01 <- list()
index1 <- list()
lambda1 <- 0
knn1 <- 0
beta1 <- list()
IC1 <- rep(Inf, length=num.wn)
vn1<-numeric(num.wn)
index2 <- list()
lambda2 <- 0
knn2 <- 0
beta2 <- list()
IC2 <- rep(Inf, length=num.wn)
vn2<-numeric(num.wn)
indexes.beta.nonnull2 <- list()
for (w in 1:num.wn) {
	message(w, "/", num.wn)
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
	step.1 <-sfpl.kNN.fit(x=x[train.1,], z=z[train.1,index.1], y=y[train.1], 
						  semimetric=semimetric, q=q, nknot=nknot,
						  lambda.min=lambda.min, lambda.min.h=lambda.min.pn.high, lambda.min.l=lambda.min.pn.low, factor.pn=factor.pn,
						  nlambda=nlambda, lambda.seq=NULL, vn=vn, nfolds=nfolds, seed=seed, 
						  knearest=knearest, min.knn=min.knn, max.knn=max.knn, step=step,
						  range.grid=range.grid, kind.of.kernel=kind.of.kernel,
						  criterion=criterion, penalty=penalty, 
						  max.iter=max.iter)
	beta <- step.1$beta.est
	index.X.pen <- step.1$indexes.beta.nonnull
	beta1[[w]] <- beta
	index1[[w]] <- index.1[index.X.pen]
	lambda1[w] <- step.1$lambda.opt 
	knn1[w] <- step.1$k.opt
	IC1[w] <- step.1$IC
	vn1[w]<-step.1$vn.opt
	if (sum(index.X.pen)==0){		
		beta2[[w]] <- NaN
		MSEP2[w] <- NaN
		indexes.beta.nonnull2[[w]] <- NaN
		IC2[w] <- NaN
		knn2[w] <- NaN
		index2[[w]] <- NaN
		next
	}		
	index.2 <- 0
	index.infl <- (1:wn[w])[index.X.pen]
	for (j in index.infl) index.2 <- c(index.2, index.wn[[j]])
	index.2 <- index.2[-1]
	index2[[w]] <- index.2
	step.2 <-sfpl.kNN.fit(x=x[train.2,], z=z[train.2,index.2], y=y[train.2], 
						  semimetric=semimetric, q=q, nknot=nknot,
						  lambda.min=lambda.min, lambda.min.h=lambda.min.pn.high, lambda.min.l=lambda.min.pn.low, factor.pn=factor.pn,
						  nlambda=nlambda, lambda.seq=NULL, vn=vn, nfolds=nfolds, seed=seed, 
						  knearest=knearest, min.knn=min.knn, max.knn=max.knn, step=step,
						  range.grid=range.grid, kind.of.kernel=kind.of.kernel,
						  criterion=criterion, penalty=penalty, 
						  max.iter=max.iter)
	beta2[[w]] <- step.2$beta.est
	indexes.beta.nonnull2[[w]] <- indexes.beta[index.2][step.2$indexes.beta.nonnull]
	lambda2[w] <- step.2$lambda.opt
	knn2[w] <- step.2$k.opt
	IC2[w] <- step.2$IC
	vn2[w]<- step.2$vn.opt			
} 
index.w.opt<- order(IC2)[1]
beta.red<-beta2[[index.w.opt]]
w.opt<-wn[index.w.opt]
IC<-IC2[index.w.opt]
vn.opt<-vn2[index.w.opt]
lambda.opt<-lambda2[index.w.opt]
k.opt<-knn2[index.w.opt]
indexes.beta.nonnull<-indexes.beta.nonnull2[[index.w.opt]]
beta.red.nonnull<-beta.red[beta.red!=0]
beta.est<-numeric(pn)
beta.est[indexes.beta.nonnull]<-beta.red.nonnull
ww<-H.fnp.kNN(x=x,k=k.opt,semimetric=semimetric,kind.of.kernel=kind.of.kernel,q=q,range.grid=range.grid,nknot=nknot)
yhp<-z%*%beta.est+ww%*%(y-z%*%beta.est)
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta.est, indexes.beta.nonnull=indexes.beta.nonnull,w.opt=w.opt,k.opt=k.opt,
		  lambda.opt=lambda.opt,IC=IC,vn.opt=vn.opt,
		  beta=beta2, IC2=IC2, lambda2=lambda2, knn2=knn2, index02=index2, indexes.beta.nonnull=indexes.beta.nonnull2, 
		  beta1=beta1, IC1=IC1, lambda1=lambda1, knn1=knn1, index01=index01, index1=index1, 
		  call=call,x=x,y=y,z=z,train.1=train.1,train.2=train.2,
		  kind.of.kernel=kind.of.kernel, range.grid=range.grid,nknot=nknot,
		  semimetric=semimetric, q=q,
		  penalty=penalty,criterion=criterion,max.iter=max.iter,wn=wn,
		  seed=seed,nlambda=nlambda,lambda.min.pn.low=lambda.min.pn.low, lambda.min.pn.high=lambda.min.pn.high, factor.pn=factor.pn,
		  nfolds=nfolds,group=group,vn=vn,
		  step=step,knearest=knearest,min.knn=min.knn,max.knn=max.knn)
class(out)<-"PVS.kNN"
out
}











