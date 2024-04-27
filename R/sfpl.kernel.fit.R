sfpl.kernel.fit<- function(x, z, y,
semimetric="deriv", q=NULL,
min.q.h=0.05, max.q.h=0.5, h.seq = NULL, num.h =10,
range.grid=NULL,  kind.of.kernel="quad", nknot=NULL, 
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, lambda.seq=NULL,vn=ncol(z), nfolds=10, seed=123,   
criterion="GCV",  
penalty="grSCAD", 
max.iter=1000)
{
if (!is.matrix(z)) z <- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
n <- nrow(z)
num.vn <- length(vn)
min.quantile.h<-min.q.h
max.quantile.h<-max.q.h
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
pn <- ncol(z)
indexes.beta <- 1:pn
p <- ncol(x)
kind.of.semimetric <- paste("semimetric.", semimetric, sep = "")
if (is.null(q)) q<-ifelse(semimetric=="deriv",0,2)
if (is.null(lambda.min)) {
	if (is.null(lambda.min.pn.high)) lambda.min.pn.high <- 0.05
	if (is.null(lambda.min.pn.low)) lambda.min.pn.low <- 1e-4
	lambda.min={if (nrow(z) > (factor.pn*ncol(z))) lambda.min.pn.low else lambda.min.pn.high} 
}
if (!(is.null(h.seq))) num.h <- length(h.seq)
if (is.null(range.grid)) range.grid <- c(1,p)
lambda2 <- 0
h2 <- 0
beta2 <- list()
index2 <- list()
lambdas<-list()
IC2 <- rep(Inf, length=num.vn)
posicion.lambda.mean.dt <- matrix(0,num.vn,2)
h.min.opt.max.mopt <- matrix(NA,num.vn,3)
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
	XX <- array(0,c(n,pn,num.h))
	aux.Yhat.hseq <- fun.kernel(y=y, x=x, pred=x, semimetric=semimetric, q=q, min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h, h.seq = h.seq, num.h = num.h, range.grid=range.grid, nknot=nknot, kind.of.kernel=kind.of.kernel)
	yy <- y- aux.Yhat.hseq$yhat
	for (j in 1:pn) XX[,j,]<- z[,j]- fun.kernel(y=z[,j], x=x, pred=x, semimetric=semimetric, q=q, min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h, h.seq = h.seq, num.h= num.h, range.grid=range.grid, nknot=nknot, kind.of.kernel=kind.of.kernel)$yhat
	h.seq <- aux.Yhat.hseq$h.seq
	lambda.s <- 0
	IC.s <- 0
	if (criterion != "k-fold-CV") {
		for (s in 1:num.h) {
			if (is.null(lambda.seq)) aux0 <-  try(grpreg(X=as.matrix(XX[,,s]), y=yy[,s], group=group, lambda.min=lambda.min, nlambda=nlambda, penalty=penalty, max.iter=max.iter), silent=TRUE)
			else aux0 <-  try(grpreg(X=as.matrix(XX[,,s]), y=yy[,s], group=group, lambda=lambda.seq, penalty=penalty, max.iter=max.iter), silent=TRUE)
			if (inherits(aux0,"try-error")) {
				lambda.s[s] <- NaN
				IC.s[s] <- NaN
				next
			}
			select<-grpreg::select
			aux <- select(obj=aux0, criterion=criterion)
			index.finite <- is.finite(aux$IC)	
			IC.finite <- aux$IC[index.finite]
			opt <- order(IC.finite)[1]
			IC.s[s] <- IC.finite[opt]
			lambda.s[s] <- aux0$lambda[index.finite][opt]
		}
		s.opt <- order(IC.s)[1]
		lambda2[v] <- lambda.s[s.opt]
		lambdas[[v]]<-lambda.s
		h2[v] <- h.seq[s.opt]
		IC2[v] <- IC.s[s.opt]
		aux <- try(grpreg(X=as.matrix(XX[,,s.opt]), y=yy[,s.opt], group=group, lambda=lambda2[v], penalty=penalty, max.iter=max.iter), silent=TRUE)
		if (inherits(aux,"try-error")) next
		beta2[[v]] <- aux$beta[-1]
		index2[[v]] <- indexes.beta[beta2[[v]]!=0]
	} 
	else {
		ss <- 0
		for (s in 1:num.h) {  
			if (is.null(lambda.seq)) aux <-  try(cv.grpreg(lambda.min=lambda.min, nlambda=nlambda, X=as.matrix(XX[,,s]), y=yy[,s], group=group, penalty=penalty, 
														   nfolds=nfolds, max.iter=max.iter, seed=seed), silent=TRUE)
			else aux <-  try(cv.grpreg(lambda=lambda.seq, X=as.matrix(XX[,,s]), y=yy[,s], group=group, penalty=penalty, nfolds=nfolds, max.iter=max.iter, seed=seed), silent=TRUE)
			if (inherits(aux,"try-error")) {
				lambda.s[s] <- NaN
				IC.s[s] <- NaN
				next
			}
			ss <- ss + 1				
			index.finite <- aux$cve >=0
			IC.finite <- aux$cve[index.finite]
			opt <- order(IC.finite)[1]
			IC.s[s] <- IC.finite[opt]
			lambda.s[s] <- aux$lambda[index.finite][opt]
			posicion.s <- order(aux$lambda[index.finite]==lambda.s[s], decreasing=TRUE)[1]
			posicion.lambda.mean.dt[v,] <- posicion.lambda.mean.dt[v,] + c(posicion.s, posicion.s^2)
		}
		posicion.lambda.mean.dt[v,1] <- posicion.lambda.mean.dt[v,1]/ss
		posicion.lambda.mean.dt[v,2] <- sqrt(posicion.lambda.mean.dt[v,2]/ss - posicion.lambda.mean.dt[v,1]^2)
		s.opt <- order(IC.s)[1]
		lambda2[v] <- lambda.s[s.opt]
		h2[v] <- h.seq[s.opt]
		IC2[v] <- IC.s[s.opt]
		aux <- try(grpreg(X=as.matrix(XX[,,s.opt]), y=yy[,s.opt], group=group, lambda=lambda2[v], penalty=penalty, max.iter=max.iter), silent=TRUE)
		if (inherits(aux,"try-error")) next
		beta2[[v]] <- aux$beta[-1]
		index2[[v]] <- indexes.beta[beta2[[v]]!=0]
	} 
	h.min.opt.max.mopt[v,] <- c(min(h.seq),h2[v], max(h.seq))
} 
ind.vn<-order(IC2)[1]
vn.opt<-vn[ind.vn]
ww<-H.fnp.kernel(x=x,h=h2[ind.vn],semimetric=semimetric,q=q,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot)
yhp<-z%*%beta2[[ind.vn]]+ww%*%(y-z%*%beta2[[ind.vn]])
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta2[[ind.vn]],indexes.beta.nonnull=index2[[ind.vn]], h.opt=h2[ind.vn],
		  lambda.opt=lambda2[ind.vn], IC=IC2[ind.vn], h.min.opt.max.mopt=h.min.opt.max.mopt[ind.vn,],  
		  vn.opt=vn.opt,lambda.opt.h=lambdas,
		  call=call,y=y,x=x,z=z, n=n, 
		  kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,max.iter=max.iter,
		  semimetric=semimetric,q=q, 
		  nfolds=nfolds,group=group,vn=vn,seed=seed, 
		  penalty=penalty, criterion=criterion,
		  num.h=num.h,h.seq=h.seq,min.quantile.h=min.quantile.h,max.quantile.h=max.quantile.h)
class(out)<-"sfpl.kernel"
return(out)
}












