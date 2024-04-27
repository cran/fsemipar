sfplsim.kernel.fit.fixedtheta<- function(x,z, y,
norm.diff,
min.quantile.h=0.05, max.quantile.h=0.5, h.seq = NULL, num.h = 10,
range.grid=NULL, kernel, 
lambda.min=NULL, lambda.min.pn.high=NULL, lambda.min.pn.low=NULL, factor.pn=1,
nlambda=100, lambda.seq=NULL,vn=ncol(z), nfolds=10, seed=123, 
criterion="GCV",  
penalty="grSCAD",  
max.iter=1000)
{
if (penalty=="grSCAD") penal <- get("SCAD")
else if (penalty=="grLasso") penal <- get("LASSOmia")
else stop("Only grLasso and grSCAD are implemented")
if (!is.matrix(z)) z <- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
num.vn <- length(vn)
n <- nrow(z)
pn <- ncol(z)
indexes.beta <- 1:pn
p <- ncol(x)
if (is.null(lambda.min)) {
	if (is.null(lambda.min.pn.high)) lambda.min.pn.high <- 0.05
	if (is.null(lambda.min.pn.low)) lambda.min.pn.low <- 1e-4
	lambda.min={if (nrow(z) > (factor.pn*ncol(z))) lambda.min.pn.low else lambda.min.pn.high} 
}
if (!(is.null(h.seq))) num.h <- length(h.seq)
lambda2 <- 0
lambda3 <- matrix(0,num.vn,3)
h2 <- 0
h3 <- matrix(0,num.vn,3)
H.SEQ <- vector("list",num.vn)
beta2 <- vector("list",num.vn)
index2 <-vector("list",num.vn)
indexes.beta.nonnull <-vector("list",num.vn)
IC2 <- rep(Inf, length=num.vn)
Q2 <- rep(Inf, length=num.vn)
MSEP2 <- rep(Inf, length=num.vn)
for (v in 1:num.vn) {
    vnv<- vn[v]
	num.veci <- trunc(pn/vnv)
	aux <- pn - vnv*num.veci
	group <- 0
	if (aux!=0) {
		for (j in 1:aux) group <- c(group, rep(j, length=num.veci+1))
		for (j in (aux+1):vnv) group <- c(group, rep(j, length=num.veci))
	}
	else for (j in 1:vnv) group <- c(group, rep(j, length=num.veci))
	group <- group[-1]
	XX <- array(0,c(n,pn,num.h))
	ajustey <- fun.kernel.fixedtheta(y=y,norm.diff=norm.diff, min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h,h.seq = h.seq, num.h = num.h,kernel=kernel)
	h.seq <- ajustey$h.seq
	H.SEQ[[v]] <- h.seq
	yy <- y-ajustey$Yhat
	for (j in 1:pn) XX[,j,]<- z[,j]-fun.kernel.fixedtheta(y=z[,j], norm.diff=norm.diff, min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h,h.seq = h.seq,num.h = num.h,
														   kernel=kernel)$Yhat
	lambda.s <- numeric(num.h)
	lambda3.s <- matrix(0,num.h,3)
	IC.s <- numeric(num.h)
	Q.s <- numeric(num.h)
	if (criterion != "k-fold-CV") {
		for (s in 1:num.h) {
			if (is.null(lambda.seq)) aux0 <-  try(grpreg(X=as.matrix(XX[,,s]), y=yy[,s], lambda.min=lambda.min, nlambda=nlambda, group=group, penalty=penalty, max.iter=max.iter), silent=TRUE)
			else aux0 <-  try(grpreg(X=as.matrix(XX[,,s]), y=yy[,s], lambda=lambda.seq, group=group, penalty=penalty, max.iter=max.iter), silent=TRUE)
			if (inherits(aux0,"try-error")) {
				lambda.s[s] <- NaN
				lambda3.s[s,] <- NaN
				IC.s[s] <- NaN
				next
			}
			select<-grpreg::select
			aux <- select(obj=aux0, criterion=criterion)    
			index.finite <- is.finite(aux$IC)			
			IC.finite <- aux$IC[index.finite]
			opt <- which.min(IC.finite)
			IC.s[s] <- IC.finite[opt]
			lambda.s[s] <- aux0$lambda[index.finite][opt]
			lambda3.s[s,] <- c(min(aux0$lambda[index.finite]), lambda.s[s], max(aux0$lambda[index.finite]))
			Q.s[s] <- 0.5*aux0$loss[index.finite][opt] 
		}
	} 
	else {
		for (s in 1:num.h) {
			if (is.null(lambda.seq)) aux <-  try(cv.grpreg(lambda.min=lambda.min, nlambda=nlambda, X=as.matrix(XX[,,s]), y=yy[,s], group=group, penalty=penalty, nfolds=nfolds, max.iter=max.iter, seed=seed), silent=TRUE)
			else aux <-  try(cv.grpreg(lambda=lambda.seq, X=as.matrix(XX[,,s]), y=yy[,s], group=group, penalty=penalty, nfolds=nfolds, max.iter=max.iter, seed=seed), silent=TRUE)
			if (inherits(aux,"try-error")) {
				lambda.s[s] <- NaN
				lambda3.s[s,] <- NaN
				IC.s[s] <- NaN
				next
			}
			index.finite <- aux$cve >=0
			IC.finite <- aux$cve[index.finite]
			opt <- which.min(IC.finite)
			IC.s[s] <- IC.finite[opt]
			lambda.s[s] <- aux$lambda[index.finite][opt]
			lambda3.s[s,] <- c(min(aux$lambda[index.finite]), lambda.s[s], max(aux$lambda[index.finite]))
			Q.s[s] <- 0.5*aux$fit$loss[index.finite][opt] 
		}		
	} 
	s.opt <- which.min(IC.s)
	lambda2[v] <- lambda.s[s.opt]
	lambda3[v,] <- lambda3.s[s.opt,]
	h2[v] <- h.seq[s.opt]
	h3[v,] <- c(min(H.SEQ[[v]]), h2[v], max(H.SEQ[[v]]))
	IC2[v] <- IC.s[s.opt]
	aux <- try(grpreg(X=as.matrix(XX[,,s.opt]), y=yy[,s.opt], group=group, lambda=lambda2[v], penalty=penalty, max.iter=max.iter), silent=FALSE)
	if (inherits(aux,"try-error")) next
	beta2[[v]] <- aux$beta
	indexes.beta.nonnull[[v]] <-indexes.beta[beta2[[v]][-1]!=0]
	beta2[[v]] <- beta2[[v]][-1] 
	sd.XX.s <- apply(XX[,,s.opt], 2, "sd")
	beta.b <-  sd.XX.s*beta2[[v]]
	Q2[v] <- Q.s[s.opt] + n*sum(penal(beta.b,lambda2[v]))
} 
ind.vn<-which.min(IC2)
vn.opt<-vn[ind.vn]
list(beta=beta2, indexes.beta.nonnull=indexes.beta.nonnull, IC=IC2[ind.vn], lambda=lambda2[ind.vn], lambda.min.opt.max=lambda3[ind.vn,],
	 h=h2[ind.vn], h.min.opt.max=h3[ind.vn,], H.SEQ=H.SEQ, Q=Q2[ind.vn], criterion=criterion, penalty=penalty,vn.opt=vn.opt,ind.vn=ind.vn)
}


















