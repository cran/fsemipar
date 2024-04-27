PVS.fit<- function(z, y, train.1=NULL, train.2=NULL,  
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, vn=ncol(z), nfolds=10, seed=123,wn=c(10,15,20),range.grid=NULL,
criterion="GCV", 
penalty="grSCAD", 
max.iter=1000)
{
if (!is.matrix(z)) z<- as.matrix(z)
onez <- cbind(rep(1, nrow(z)), z)
n<-nrow(z)
if(is.null(train.1)) train.1<-1:ceiling(n/2)
if(is.null(train.2)) train.2<-(ceiling(n/2)+1):n
n1 <- length(train.1)
n2 <- length(train.2)
num.wn <- length(wn)
pn <- ncol(z)
if(is.null(range.grid)) range.grid=c(1,pn)
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
indexes.beta <- 1:pn
index01 <- list()
index1 <- list()
lambda1 <- 0
beta1 <- list()
IC1 <- rep(Inf, length=num.wn)
vn1<-numeric(num.wn)
index2 <- list()
lambda2 <- 0
beta2 <- list()
IC2 <- rep(Inf, length=num.wn)
vn2<-numeric(num.wn)
indexes.beta.nonnull2 <- list()
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
	step.1 <- lm.pels.fit(z=z[train.1,index.1], y=y[train.1], 
						  lambda.min=lambda.min, lambda.min.h=lambda.min.pn.high, 
						  lambda.min.l=lambda.min.pn.low, factor.pn=factor.pn,
						  nlambda=nlambda, lambda.seq=NULL, vn=vn, nfolds=nfolds, seed=seed, 
						  criterion=criterion, penalty=penalty, 
						  max.iter=max.iter)
	beta <- step.1$beta.est
	index.X.pen <- step.1$indexes.beta.nonnull
	beta1[[w]] <- beta
	index1[[w]] <- index.1[index.X.pen]
	lambda1[w] <- step.1$lambda.opt 
	IC1[w] <- step.1$IC
	vn1[w]<- step.1$vn.opt
	if (sum(index.X.pen)==0){		
		beta2[[w]] <- NaN
		MSEP2[w] <- NaN 
		indexes.beta.nonnull2[[w]] <- NaN
		IC2[w] <- NaN 
		index2[[w]] <- NaN
		next
	}
	index.2 <- 0
	index.infl <- (1:wn[w])[index.X.pen]
	for (j in index.infl) index.2 <- c(index.2, index.wn[[j]])
	index.2 <- index.2[-1]
	index2[[w]] <- index.2
	step.2 <- lm.pels.fit(z=z[train.2,index.2], y=y[train.2], 
						  lambda.min=lambda.min, lambda.min.h=lambda.min.pn.high, 
						  lambda.min.l=lambda.min.pn.low, factor.pn=factor.pn,
						  nlambda=nlambda, lambda.seq=NULL, vn=vn, nfolds=nfolds, seed=seed, 
						  criterion=criterion, penalty=penalty, 
						  max.iter=max.iter)
	beta2[[w]] <- step.2$beta.est
	indexes.beta.nonnull2[[w]] <- indexes.beta[index.2][step.2$indexes.beta.nonnull]
	lambda2[w] <- step.2$lambda.opt
	IC2[w] <- step.2$IC
	vn2[w]<- step.2$vn.opt
} 
index.w.opt<- order(IC2)[1]
beta.red<-beta2[[index.w.opt]]
w.opt<-wn[index.w.opt]
IC<-IC2[index.w.opt]
vn.opt<-vn2[index.w.opt]
lambda.opt<-lambda2[index.w.opt]
indexes.beta.nonnull<-indexes.beta.nonnull2[[index.w.opt]]
beta.red.nonnull<-beta.red[beta.red!=0]
beta.est<-numeric(pn+1)
beta.est[1]<-beta.red.nonnull[1]
beta.est[indexes.beta.nonnull]<-beta.red.nonnull[2:length(beta.red.nonnull)]
yhp<-as.matrix(onez)%*%beta.est
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta.est,indexes.beta.nonnull=indexes.beta.nonnull2, IC=IC, lambda.opt=lambda.opt,vn.opt=vn.opt,
		  beta2=beta2,  IC2=IC2,   lambda2=lambda2, indexes.beta.nonnull=indexes.beta.nonnull2,  index02=index2,
		  beta1=beta1, IC1=IC1, lambda1=lambda1, index01=index01, index1=index1, 
		  call=call,y=y,z=z,train.1=train.1, train.2=train.2,range.grid=range.grid,
		  penalty=penalty, criterion=criterion,max.iter=max.iter,wn=wn)
class(out)<-"PVS"
out
}