sfplsim.kNN.fit.fixedtheta<- function(x, z, y, norm.diff,  
knearest=NULL, min.knn=2, max.knn=NULL,step=NULL, kernel,
lambda.min=NULL, lambda.min.pn.high=NULL, lambda.min.pn.low=NULL, factor.pn=1,
nlambda=100, lambda.seq=lambda.seq,vn=ncol(z), nfolds=10, seed=123, 
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
if (is.null(max.knn)) max.knn <- n%/%2
if (is.null(knearest)) {
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
num.knn <- length(knearest)
lambda2 <- 0
lambda3 <- matrix(0,num.vn,3)
knn2 <- 0
knn3 <- matrix(0,num.vn,3)
beta2 <- list()
index2 <- list()
indexes.beta.nonnull <- list()
IC2 <- rep(Inf, length=num.vn)
Q2 <- rep(Inf, length=num.vn)
MSEP2 <- rep(Inf, length=num.vn)
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
# GROUP
	XX<- array(0,c(n,pn,num.knn))
	yy <- y-fun.kNN.fixedtheta(y=y, norm.diff=norm.diff, knearest=knearest, kernel=kernel)$yhat
	for (j in 1:pn) XX[,j,]<- z[,j]- fun.kNN.fixedtheta(y=z[,j],norm.diff=norm.diff, knearest=knearest, kernel=kernel)$yhat
	lambda.s <- 0
	lambda3.s <- matrix(0,num.knn,3)
	IC.s <- 0
	Q.s <- 0
	if (criterion != "k-fold-CV") {
		for (s in 1:num.knn) {
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
            index.finite<-is.finite(aux$IC)
			IC.finite <- aux$IC[index.finite]
			opt <- order(IC.finite)[1]
			IC.s[s] <- IC.finite[opt]
			lambda.s[s] <- aux0$lambda[index.finite][opt]
			lambda3.s[s,] <- c(min(aux0$lambda[index.finite]), lambda.s[s], max(aux0$lambda[index.finite]))
			Q.s[s] <- 0.5*aux0$loss[index.finite][opt] 
		}
	} 
	else {
		for (s in 1:num.knn) {
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
			opt <- order(IC.finite)[1]
			IC.s[s] <- IC.finite[opt]
			lambda.s[s] <- aux$lambda[index.finite][opt]
			lambda3.s[s,] <- c(min(aux$lambda[index.finite]), lambda.s[s], max(aux$lambda[index.finite]))
			Q.s[s] <- 0.5*aux$fit$loss[index.finite][opt] 
		} 		
	} 
	s.opt <- order(IC.s)[1]
	lambda2[v] <- lambda.s[s.opt]
	lambda3[v,] <- lambda3.s[s.opt,]
	knn2[v] <- knearest[s.opt]
	knn3[v,] <- c(min(knearest), knn2[v], max(knearest))
	IC2[v] <- IC.s[s.opt]
	aux <- try(grpreg(X=as.matrix(XX[,,s.opt]), y=yy[,s.opt], group=group, lambda=lambda2[v], penalty=penalty, max.iter=max.iter), silent=TRUE)	
	if (inherits(aux,"try-error")) next
	beta2[[v]] <- aux$beta
	indexes.beta.nonnull[[v]] <-indexes.beta[beta2[[v]][-1]!=0]
	beta2[[v]] <- beta2[[v]][-1] 
	sd.XX.s <- apply(XX[,,s.opt], 2, "sd")
	beta.b <-  sd.XX.s*beta2[[v]]
	Q2[v] <- Q.s[s.opt] + n*sum(penal(beta.b,lambda2[v]))
} 
ind.vn<-order(IC2)[1]
vn.opt<-vn[ind.vn]
list(beta=beta2, indexes.beta.nonnull=indexes.beta.nonnull, IC=IC2[ind.vn], lambda=lambda2[ind.vn], lambda.min.opt.max=lambda3[ind.vn,], 
     knn.min.opt.max=knn3[ind.vn,], knn=knn2[ind.vn], knearest=knearest, Q=Q2[ind.vn], criterion=criterion, penalty=penalty,vn.opt=vn.opt,ind.vn=ind.vn)
}












