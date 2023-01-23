sfplsim.kernel.fit<- function(x,z,y, 
seed.coeff=c(-1,0,1), order.Bspline=3, nknot.theta=3,   t0=NULL, 
min.q.h=0.05, max.q.h=0.5, h.seq = NULL, num.h = 10,
range.grid=NULL, kind.of.kernel="quad",nknot=NULL,
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, lambda.seq=NULL,vn=ncol(z), nfolds=10, seed=123,
criterion=c("GCV", "BIC", "AIC", "k-fold-CV"),   
penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"), 
max.iter=1000)
{
if (!is.matrix(z)) z <- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
n <- nrow(z)
pn <- ncol(z)
indexes.beta <- 1:pn
min.quantile.h<-min.q.h
max.quantile.h<-max.q.h
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
p <- ncol(x)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
if (is.null(lambda.min)) {
	if (is.null(lambda.min.pn.high)) lambda.min.pn.high <- 0.05
	if (is.null(lambda.min.pn.low)) lambda.min.pn.low <- 1e-4
	lambda.min={if (nrow(z) > (factor.pn*ncol(z))) lambda.min.pn.low else lambda.min.pn.high} 
}
if (!(is.null(h.seq))) num.h <- length(h.seq)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(t0)) t0 <- mean(range.grid)
dim.base.theta <- order.Bspline + nknot.theta
THETA.seq <- permutations(length(seed.coeff), dim.base.theta, seed.coeff, repeats.allowed=TRUE)
THETA.seq <- THETA.seq[(apply(abs(THETA.seq) , 1,sum) != 0)  & (!is.na(apply(abs(THETA.seq), 1,sum))), ]
THETA.seq.normalizado <- normaliza(coef=THETA.seq, range.grid=range.grid, t0=t0, order.Bspline=order.Bspline, nknot.theta=nknot.theta)
num.norm <- nrow(THETA.seq.normalizado)
Q <- 0
h <- 0
lambda <- 0
lambda3 <- matrix(0,num.norm,3)
h3 <- matrix(0,num.norm,3)
IC <- 0
H.SEQ <- list()
h.opt.2 <- NA
h2.min.opt.max.mopt <- rep(NA,3)
H2.SEQ <- NA
prueba<-list()
for(m in 1:num.norm) {
	message(m,"/",num.norm)
	aux <- sfplsim.kernel.fit.fixedtheta(x=x, z=z, y=y, theta=THETA.seq.normalizado[m,], 
										 order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta,
										 lambda.min=lambda.min, nlambda=nlambda, lambda.seq=lambda.seq, vn=vn, nfolds=nfolds, seed=seed, 
										 min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h, h.seq = h.seq, num.h = num.h,
										 range.grid=range.grid, kind.of.kernel=kind.of.kernel, criterion=criterion, penalty=penalty, max.iter=max.iter)
	Q[m] <- aux$Q
	IC[m] <- aux$IC
	h[m] <- aux$h
	lambda[m] <- aux$lambda
	lambda3[m,] <- aux$lambda.min.opt.max
	h3[m,] <- aux$h.min.opt.max
	H.SEQ[[m]] <- aux$H.SEQ
	prueba[[m]]<-aux$prueba
}
m.opt <- order(Q)[1]
Q.opt <- Q[m.opt]
IC.opt <- IC[m.opt]
theta.opt <- THETA.seq.normalizado[m.opt,]
h.opt <- h[m.opt]
lambda.opt <- lambda[m.opt]
aux <- sfplsim.kernel.fit.fixedtheta(x=x,y=y,z=z, theta=theta.opt, h.seq=h.opt, 
									 order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta, 
									 lambda.min=lambda.min, nlambda=nlambda, lambda.seq=lambda.seq, vn=vn, nfolds=nfolds, seed=seed, 
									 range.grid=range.grid, kind.of.kernel=kind.of.kernel, criterion=criterion, penalty=penalty, max.iter=max.iter)
ww<-H.fsim.kernel(x=x,h=h.opt,theta=theta.opt,nknot.theta=nknot.theta,order.Bspline=order.Bspline,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot)
yhp<-z%*%aux$beta[[aux$ind.vn]]+ww%*%(y-z%*%aux$beta[[aux$ind.vn]])
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=aux$beta[[aux$ind.vn]], theta.est=theta.opt, indexes.beta.nonnull=aux$indexes.beta.nonnull[[aux$ind.vn]], h.opt=h.opt,
		  lambda.opt=aux$lambda,IC=aux$IC, Q.opt=aux$Q, Q=Q,
		  m.opt=m.opt, lambda.min.opt.max.mopt=lambda3[m.opt,], lambda.min.opt.max.m=lambda3, 
		  h.min.opt.max.mopt=h3[m.opt,], h.min.opt.max.m=h3, h.seq.opt=H.SEQ[[m.opt]], vn.opt=aux$vn.opt,
		  theta.seq.norm=THETA.seq.normalizado, 
		  call=call,y=y,x=x,z=z,n=n,
		  kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot, max.iter=max.iter,
		  order.Bspline=order.Bspline,nknot.theta=nknot.theta,seed.coeff=seed.coeff,t0=t0,
		  penalty=penalty, criterion=criterion,
		  num.h=num.h,h.seq=h.seq,min.quantile.h=min.quantile.h,max.quantile.h=max.quantile.h)
class(out)<-"sfplsim.kernel"
out
}







