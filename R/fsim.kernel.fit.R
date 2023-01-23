fsim.kernel.fit<- function(x,y,
seed.coeff=c(-1,0,1),  nknot.theta=3,order.Bspline=3, t0=NULL, 
min.q.h=0.05, max.q.h=0.5, h.seq=NULL, num.h=10, 
kind.of.kernel="quad", range.grid=NULL, nknot=NULL)
{
if (!is.matrix(x))  stop("x must contain a matrix")
if (!is.matrix(y)) y <- as.matrix(y)
kernel <- get(kind.of.kernel)
n <- nrow(x)
min.quantile.h<-min.q.h
max.quantile.h<-max.q.h
if (!(is.null(h.seq))) num.h <- length(h.seq)
p <- ncol(x)
length.curve.y<-ncol(y)
estimated.Y <- list()
yhat.cv <- matrix(0,n,length.curve.y)
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(t0)) t0 <- mean(range.grid)
if (is.null(nknot))  nknot <- (p - 0 - order.Bspline - 1)%/%2 
dim.base.theta <- order.Bspline + nknot.theta
THETA.seq <- permutations(length(seed.coeff), dim.base.theta, seed.coeff, repeats.allowed=TRUE)
THETA.seq <- THETA.seq[(apply(abs(THETA.seq) , 1,sum) != 0)  & (!is.na(apply(abs(THETA.seq), 1,sum))), ]
THETA.seq.normalizado <- normaliza(coef=THETA.seq, range.grid=range.grid, t0=t0, order.Bspline=order.Bspline, nknot.theta=nknot.theta)
num.norm <- nrow(THETA.seq.normalizado)
h.seq.app <- list()
h.vec <- rep(0,num.norm)									
CV.app <- 0     
for(m in 1:num.norm) {
	message(m,"/",num.norm)
	h.seq.m <- h.seq
	num.h.m <- num.h
	if (is.null(h.seq)){ 
		norm.diff.0 <- semimetric.projec(data1=x, data2=x, theta=THETA.seq.normalizado[m,], range.grid=range.grid, order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)
		norm.diff.00 <- norm.diff.0[row(norm.diff.0) > col(norm.diff.0)]
		h.seq.m <- quantile(norm.diff.00, seq(min.quantile.h, max.quantile.h, length = num.h))
		h.seq.m <- h.seq.m[h.seq.m>0]
		num.h.m <- length(h.seq.m)
		cv.hseq <- rep(0,num.h.m) 
	}
	h.seq.app[[m]] <- h.seq.m
	for(i in 1:n){       
		for (j in 1:num.h.m) {
			h <- h.seq.m[j]
			yhat1 <- fsim.kernel.test(y=y[-i],x=x[-i,], y.test=y[i],x.test=x[i,], theta=THETA.seq.normalizado[m,], h=h,  kind.of.kernel=kind.of.kernel, range.grid=range.grid, 
									  order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)$y.estimated.test
			if (length.curve.y==1) resid.hseq.2 <- (yhat1- y[i])^2
			else resid.hseq.2 <- sum((yhat1- y[i,])^2)            
			cv.hseq[j] <- cv.hseq[j] + resid.hseq.2
		} 		                 
	} 
	cv.hseq <- cv.hseq/n
	index <- order(cv.hseq)[1]
	h.opt.m <- h.seq.m[index]
	h.vec[m] <- h.opt.m
	CV.app[m] <- cv.hseq[index]
	for (i in 1:n) 
		yhat.cv[i,] <- fsim.kernel.test(y=y[-i], x=x[-i,],y.test=y[i],x.test=x[i,], theta=THETA.seq.normalizado[m,], h=h.opt.m,  kind.of.kernel=kind.of.kernel, range.grid=range.grid,
										order.Bspline=order.Bspline, nknot=nknot, nknot.theta=nknot.theta)$y.estimated.test
		estimated.Y[[m]] <- yhat.cv
} 
m.opt <- order(CV.app)[1]
theta.opt <- THETA.seq.normalizado[m.opt,]
h.opt <- h.vec[m.opt]
h.seq.opt <- h.seq.app[[m.opt]]
H<-H.fsim.kernel(x=x,h=h.opt,theta=theta.opt,nknot.theta=nknot.theta,order.Bspline=order.Bspline,kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot)
CV.opt<-CV.app[m.opt]
fitted.values<-H%*%y
res<-y-drop(fitted.values)
df<-sum(diag(H))
sigma.res<-sum(res^2)/(n-df)
r.2<-1-sum(res^2)/sum((y-mean(y))^2)
call<-match.call()
out<-list(fitted.values=fitted.values, residuals=res,theta.est=theta.opt, h.opt=h.opt,r.squared=r.2,var.res=sigma.res,df=df,
		  yhat.cv=estimated.Y[[order(CV.app)[1]]],CV.opt=CV.opt,CV.values=CV.app,H=H,
		  m.opt=m.opt, theta.seq.norm=THETA.seq.normalizado,  h.seq=h.seq.opt, 
		  call=call,y=y,x=x,n=n,
		  kind.of.kernel=kind.of.kernel,range.grid=range.grid,nknot=nknot,
		  order.Bspline=order.Bspline,nknot.theta=nknot.theta)
class(out)<-"fsim.kernel"
out
} 



