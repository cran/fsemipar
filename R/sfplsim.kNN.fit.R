sfplsim.kNN.fit<- function(x,z,y, 
seed.coeff=c(-1,0,1), order.Bspline=3,  nknot.theta=3, 
knearest=NULL, min.knn=2, max.knn=NULL, step=NULL,
range.grid=NULL,kind.of.kernel="quad", nknot=NULL,
lambda.min=NULL, lambda.min.h=NULL, lambda.min.l=NULL, factor.pn=1,
nlambda=100, lambda.seq=NULL, vn=ncol(z), nfolds=10, seed=123, 
criterion="GCV", penalty="grSCAD", max.iter=1000,n.core=NULL)
{
if (!is.matrix(z)) z <- as.matrix(z)
if (!is.matrix(x)) x <- as.matrix(x)
n<-nrow(z)
pn <- ncol(z)
lambda.min.pn.high<-lambda.min.h
lambda.min.pn.low<-lambda.min.l
indexes.beta <- 1:pn
p <- ncol(x)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
if (is.null(lambda.min)) {
	if (is.null(lambda.min.pn.high)) lambda.min.pn.high <- 0.05
	if (is.null(lambda.min.pn.low)) lambda.min.pn.low <- 1e-4
	lambda.min={if (nrow(z) > (factor.pn*ncol(z))) lambda.min.pn.low else lambda.min.pn.high} 
}
if (is.null(max.knn)) max.knn <- n%/%5
if (is.null(knearest)) {
	if (is.null(step)) step <- ceiling(n/100)
	if(step == 0) step <- 1
	knearest <- seq(from =min.knn, to = max.knn, by = step)
}
num.knn <- length(knearest)
if (is.null(range.grid)) range.grid <- c(1,p)
t0 <- mean(range.grid)
kernel<-get(kind.of.kernel)
dim.base.theta <- order.Bspline + nknot.theta
THETA.seq <- permutations(length(seed.coeff), dim.base.theta, seed.coeff, repeats.allowed=TRUE)
THETA.seq <- THETA.seq[(rowSums(abs(THETA.seq)) != 0)  & (!is.na(rowSums(abs(THETA.seq)))), ]
a <- range.grid[1]
b <- range.grid[2]
Knot.theta<-seq(a, b, length = nknot.theta + 2)[ - c(1, nknot.theta + 2)]
delta.theta<-sort(c(rep(c(a, b),order.Bspline), Knot.theta))     
point.gauss <- c(-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142)
weight.gauss <- c(0.1713244924, 0.360761573, 0.4679139346, 0.4679139346,0.360761573, 0.1713244924)
x.gauss<- 0.5 * ((b + a) + (b - a) * point.gauss)
lx.gauss<- length(x.gauss)
Bspline.g.theta<-splineDesign(delta.theta, x.gauss, order.Bspline)
H.theta <-t(Bspline.g.theta)%*%(Bspline.g.theta*(weight.gauss*0.5*(b-a)))
x0 <- seq(a, b, length = p)
Knot <-seq(a, b, length = nknot + 2)[ - c(1, nknot + 2)]  
delta <-sort(c(rep(c(a, b),order.Bspline), Knot))
Bspline <-splineDesign(delta,x0,order.Bspline)
Cmat <-crossprod(Bspline)
Dmat1 <-crossprod(Bspline, t(x))
coef.mat1 <-symsolve(Cmat, Dmat1)
Bspline.theta<-splineDesign(delta.theta,x0,order.Bspline)
Bspline.g<-splineDesign(delta, x.gauss, order.Bspline)
H.x <-t(Bspline.g)%*%(Bspline.g*(weight.gauss*0.5*(b-a)))
dim<-ncol(THETA.seq)
per<-nrow(THETA.seq)
theta.normal<-matrix(0,nrow=per,ncol=dim)
for(i in 1:per){
	coefi=THETA.seq[i,]
	prod.coef <-outer(coefi,coefi, "*")
	norm=sqrt(sum(H.theta*prod.coef))
	theta.normal[i,]=coefi/norm
}
Bspline.t0.theta<-splineDesign(delta.theta, t0, order.Bspline)
theta.t0<-theta.normal%*%t(Bspline.t0.theta)
pos<-which(theta.t0>0) 
THETA.seq.normalizado<-theta.normal[pos,]		
num.norm <- nrow(THETA.seq.normalizado)
Q <- 0
k <- 0
lambda <- 0
lambda3 <- matrix(0,num.norm,3)
knn3 <- matrix(0,num.norm,3)
IC <- 0
beta<-list()
indexes.beta.nonnull<-list()
vn.m<-numeric(num.norm)
ind<-numeric(num.norm)
if(is.null(n.core)) n.core<-availableCores(omit=1)
    cl <- makeCluster(n.core)
    registerDoParallel(cl)
	m<-NULL
	results <- foreach(m = 1:num.norm,.packages=c("grpreg")) %dopar% { 
	theta<-THETA.seq.normalizado[m,]
	theta.rec<-crossprod(t(Bspline.theta),theta) 
	Dmat.theta<-crossprod(Bspline,theta.rec)
	Theta.coef<-symsolve(Cmat,Dmat.theta)
	theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef})
	coef <-t(H.x%*%theta.x1)
	projec1<-rowSums(coef)
	semimetric <- outer(projec1,projec1,"-")
	norm.diff.0<-abs(semimetric)
	aux <- sfplsim.kNN.fit.fixedtheta(x=x, z=z, y=y, norm.diff=norm.diff.0, 
						step=step,lambda.min=lambda.min, nlambda=nlambda, lambda.seq=lambda.seq, vn=vn, nfolds=nfolds, seed=seed, knearest=knearest, min.knn=min.knn, max.knn=max.knn, kernel=kernel,
					criterion=criterion, penalty=penalty, max.iter=max.iter)
	Q[m] <- aux$Q
	IC[m] <- aux$IC
	k[m] <- aux$knn
	lambda[m] <- aux$lambda
	lambda3[m,] <- aux$lambda.min.opt.max
	knn3[m,] <- aux$knn.min.opt.max
	beta[[m]]<-aux$beta
	indexes.beta.nonnull[[m]]=aux$indexes.beta.nonnull
	ind[m]<-aux$ind.vn
	vn.m[m]<-aux$vn.opt
    return(list(Q = Q[m], 
	IC = IC[m], k = k[m], lambda = lambda[m],
    lambda3 = lambda3[m,], knn3 = knn3[m,], beta = beta[[m]],
    indexes.beta.nonnull = indexes.beta.nonnull[[m]], ind = ind[m],
    vn.m = vn.m[m]))	
}
stopCluster(cl)
Q <- sapply(results, function(r) r$Q)
m.opt <- which.min(Q)
Q.opt <- Q[m.opt]
IC.opt <- results[[m.opt]]$IC 
lambda.opt <- results[[m.opt]]$lambda
lambda3 <- results[[m.opt]]$lambda3
k.opt <- results[[m.opt]]$k
knn3 <- results[[m.opt]]$knn3
ind <- results[[m.opt]]$ind
vn.opt <- results[[m.opt]]$vn.m
beta.opt <- results[[m.opt]]$beta[[ind]]
indexes.beta.nonnull <- results[[m.opt]]$indexes.beta.nonnull[[ind]]
theta.opt <- THETA.seq.normalizado[m.opt,]
theta.rec<-crossprod(t(Bspline.theta),theta.opt) 
Dmat.theta<-crossprod(Bspline,theta.rec)
Theta.coef<-symsolve(Cmat,Dmat.theta)
theta.x1<-sapply(1:ncol(coef.mat1),function(i){coef.mat1[,i]*Theta.coef})
coef <-t(H.x%*%theta.x1)
projec1<-rowSums(coef)
semimetric <- outer(projec1,projec1,"-")
norm.diff.1<-abs(semimetric)
H<-matrix(0,nrow=n,ncol=n)
for(i in 1:n){
	norm.diff.1i<-norm.diff.1[i,]
	norm.order <- order(norm.diff.1i)
	ind.curves2 <- norm.order[1:(k.opt + 1)]
	h <- sum(abs(norm.diff.1i[ind.curves2[k.opt:(k.opt + 1)]]))/2
	res.kernel <- kernel(norm.diff.1i[ind.curves2[ - (k.opt + 1)]]/h)
	sum.res.kernel <- sum(res.kernel)
	for(j in 1:k.opt){
		cur<-ind.curves2[j]
		if(sum.res.kernel > 0) {
			H[i,cur] <-res.kernel[j]/sum.res.kernel
		}                                                
	else H[i,cur] <- 1
	}
}		
yhp<-z%*%beta.opt+H%*%(y-z%*%beta.opt)
res<-y-yhp
call<-match.call()
out<-list(fitted.values=yhp,residuals=res,beta.est=beta.opt, theta.est=theta.opt,indexes.beta.nonnull=indexes.beta.nonnull, k.opt=k.opt, 
		  lambda.opt=lambda.opt,IC=IC.opt, Q.opt=Q.opt, Q=Q, 
		  m.opt=m.opt, lambda.min.opt.max.mopt=lambda3, 
		  knn.min.opt.max.mopt=knn3, knn.min.opt.max.m=knn3,
		  theta.seq.norm=THETA.seq.normalizado, vn.opt=vn.opt, norm.diff=norm.diff.1,ww=H,
		  call=call,y=y,x=x,z=z,n=n, 
		  kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot, max.iter=max.iter,
		  order.Bspline=order.Bspline,nknot.theta=nknot.theta,seed.coeff=seed.coeff,t0=t0,
		  penalty=penalty, criterion=criterion,
		  step=step,min.knn=min.knn, max.knn=max.knn)
class(out)<-"sfplsim.kNN"
out
}
