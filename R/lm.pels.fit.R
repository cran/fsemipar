lm.pels.fit<-function(z, y,
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, lambda.seq=NULL,vn=ncol(z), nfolds=10, seed=123,  
criterion="GCV",
penalty="grSCAD", 
max.iter=1000)
{
if (!is.matrix(z)) z <- as.matrix(z)
n <- nrow(z)
num.vn <- length(vn)
pn <- ncol(z)
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
indexes.beta <- 1:pn
onez<- cbind(rep(1, nrow(z)), z)
if (is.null(factor.pn)) factor.pn <- 1
if (is.null(lambda.min)) {
	if (is.null(lambda.min.pn.high)) lambda.min.pn.high <- 0.05
	if (is.null(lambda.min.pn.low)) lambda.min.pn.low <- 1e-4
	lambda.min={if (nrow(z) > (factor.pn*ncol(z))) lambda.min.pn.low else lambda.min.pn.high} 
}
if (is.null(nlambda)) nlambda <- 100
lambda2 <- 0
beta2 <- list()
index2 <- list()
IC2 <- rep(Inf, length=num.vn)
MSEP <- rep(Inf, length=num.vn)
for (v in 1:num.vn) {
	num.veci <- trunc(pn/vn[v])
	aux <- pn - vn[v]*num.veci	
	group <- 0
	if (aux!=0) {
		for (j in 1:aux) group <- c(group, rep(j, length=num.veci+1))
		for (j in (aux+1):vn[v]) group <- c(group, rep(j, length=num.veci))
	}
	else for (j in 1:vn[v]) group <- c(group, rep(j, length=num.veci))
	group <- group[-1]	
	if (criterion != "k-fold-CV") {
		if (is.null(lambda.seq)) aux0 <-  try(grpreg(X=as.matrix(z), y=y, group=group, lambda.min=lambda.min, nlambda=nlambda, penalty=penalty, max.iter=max.iter), silent=TRUE)
		else aux0 <-  try(grpreg(X=as.matrix(z), y=y, group=group, lambda=lambda.seq, penalty=penalty, max.iter=max.iter), silent=TRUE)
		if (inherits(aux0,"try-error")) stop("the application of the selector fails over train")
		select<-grpreg::select
		aux <- select(obj=aux0, criterion=criterion)
		index.finite <- is.finite(aux$IC)
		IC.finite <- aux$IC[index.finite]
		opt <- order(IC.finite)[1]
		IC <- IC.finite[opt]
		lambda <- aux0$lambda[index.finite][opt]					
		lambda2[v] <- lambda
		IC2[v] <- IC
		aux <- try(grpreg(X=as.matrix(z), y=y, group=group, lambda=lambda2[v], penalty=penalty, max.iter=max.iter), silent=TRUE)
		if (inherits(aux,"try-error")) stop("the application of the selector with lambda.opt fails")
		beta2[[v]] <- aux$beta
		index2[[v]] <- indexes.beta[beta2[[v]]!=0]
		index2[[v]] <- indexes.beta[beta2[[v]][-1]!=0]
	} 
	else{
		if (is.null(lambda.seq)) aux <-  try(cv.grpreg(lambda.min=lambda.min, nlambda=nlambda, X=as.matrix(z), y=y, group=group, penalty=penalty, nfolds=nfolds, max.iter=max.iter, seed=seed), silent=TRUE)
		else aux <-  try(cv.grpreg(lambda=lambda.seq, X=as.matrix(z), y=y, group=group, penalty=penalty, nfolds=nfolds, max.iter=max.iter, seed=seed), silent=TRUE)
		if (inherits(aux,"try-error")) stop("the application of the selector with k-fold-CV fails")				
		index.finite <- aux$cve >=0
		IC.finite <- aux$cve[index.finite]
		opt <- order(IC.finite)[1]
		IC <- IC.finite[opt]
		lambda <- aux$lambda[index.finite][opt]
		lambda2[v] <- lambda
		IC2[v] <- IC
		aux <- try(grpreg(X=as.matrix(z), y=y, group=group, lambda=lambda2[v], penalty=penalty, max.iter=max.iter), silent=TRUE)
		if (inherits(aux,"try-error")) stop("the application of the selector with k-fold-CV and lambda.opt fails")
		index2[[v]] <- indexes.beta[beta2[[v]]!=0]
		index2[[v]] <- indexes.beta[beta2[[v]][-1]!=0]
	} 
}
ind.vn<-order(IC2)[1]
vn.opt<-vn[ind.vn] 
yhp<-onez%*%beta2[[ind.vn]]
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta2[[ind.vn]], indexes.beta.nonnull=index2[[ind.vn]], lambda.opt=lambda2[ind.vn],IC=IC2[ind.vn],vn.opt=vn.opt,
			call=call,penalty=penalty,criterion=criterion,y=y,z=z)
class(out)<-"lm.pels"
out
}



