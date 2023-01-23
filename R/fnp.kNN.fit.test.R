fnp.kNN.fit.test<- function(x,y,x.test,y.test=NULL, 
semimetric="pca", start.order.deriv.o.pca=NULL, end.order.deriv.o.pca=NULL, min.leng.interv=NULL, max.leng.interv=NULL, 
knearest=NULL, min.knn=NULL, max.knn=NULL, step=NULL, 
kind.of.kernel=NULL,range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
parameters.function <- list(x,y,x.test,y.test,semimetric=semimetric, kind.of.kernel=kind.of.kernel, start.order.deriv.o.pca=start.order.deriv.o.pca, end.order.deriv.o.pca=end.order.deriv.o.pca, min.leng.interv=min.leng.interv, max.leng.interv=max.leng.interv, range.grid=range.grid, nknot=nknot)                                                                                                                                        
p <- ncol(x)
if (is.null(kind.of.kernel)) kind.of.kernel <- "quad"
if (p==1) stop("A functional variable is needed (i.e., ncol(x)>2 is required")
if ((semimetric=="deriv") | (semimetric=="pca")){		       
	if ( (semimetric=="deriv") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 0
	if ( (semimetric=="deriv") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- 3
	if ((semimetric=="pca") & (is.null(start.order.deriv.o.pca)) )  start.order.deriv.o.pca <- 1
	if ((semimetric=="pca") & (is.null(end.order.deriv.o.pca)) )  end.order.deriv.o.pca <- p
	start.order <- start.order.deriv.o.pca
	end.order <- end.order.deriv.o.pca
}		
else if (semimetric=="interv"){
	if (is.null(min.leng.interv)) min.leng.interv <- 1
	if (is.null(max.leng.interv)) max.leng.interv <- p-1
	start.order <- 1
	end.order <- (max.leng.interv-min.leng.interv +1)*(2*p-max.leng.interv-min.leng.interv)/2
}
else stop("semimetric must be equal to: deriv, interv or pca")
if (is.null(range.grid)) range.grid <- c(1,p)
if (is.null(nknot)) {
	if (semimetric=="deriv") nknot <- (p - end.order.deriv.o.pca - 3 - 1)%/%2 
	else if (semimetric=="interv") nknot <- (p - 0 - 3 - 1)%/%2 
}
num.norm <- end.order - start.order +1
CV.app <- 1:num.norm
k.app <- 1:num.norm
for(m in start.order:end.order) {
	a <- p-min.leng.interv
	if (semimetric=="interv") {
	for (i in 1:(max.leng.interv-min.leng.interv+1)) 
		if ( ((i-1)*(2*a+2-i)/2 + 1 <= m) & (m <= (i*(2*a+1-i)/2)) ) {		                  
			ext.inf <- m - (i-1)*(2*a+2-i)/2
			ext.sup <- ext.inf + min.leng.interv + i -1
		}                                
}
	if (semimetric=="deriv") step.1 <- fnp.kNN.GCV(y=y, x=x, pred=x, q=m, knearest=knearest, min.knn=min.knn, max.knn=max.knn, step=step,
											range.grid=range.grid, nknot=nknot, kind.of.kernel=kind.of.kernel, semimetric=semimetric)
	else if (semimetric=="interv") step.1 <- fnp.kNN.GCV(y=y, x=x, pred=x, interv=c(ext.inf,ext.sup), knearest=knearest, min.knn=min.knn, max.knn=max.knn, step=step, 
											range.grid=range.grid, nknot=nknot, kind.of.kernel=kind.of.kernel, semimetric=semimetric)
	else if (semimetric=="pca") step.1 <- fnp.kNN.GCV(y=y, x=x, pred=x, q=m, knearest=knearest, min.knn=min.knn, max.knn=max.knn, step=step, 
											kind.of.kernel=kind.of.kernel, semimetric=semimetric)		
	CV.app[m + 1 - start.order] <- step.1$Mse
	k.app[m + 1 - start.order] <- step.1$knearest.opt
} # for m
m.opt <- order(CV.app)[1] +  start.order - 1
k.opt <- k.app[order(CV.app)[1]]
if (semimetric=="interv") {		
	a <- p-min.leng.interv   	
	for (i in 1:(max.leng.interv-min.leng.interv+1)) 
		if ( ((i-1)*(2*a+2-i)/2 + 1 <= m.opt) & (m.opt <= (i*(2*a+1-i)/2)) ) {
			ext.inf <- m.opt - (i-1)*(2*a+2-i)/2
			ext.sup <- ext.inf + min.leng.interv + i -1
		}
}	
if(is.vector(x.test)) x.test<- as.matrix(t(x.test))
X2.pred <- matrix(0, nrow(x.test)+1, ncol(x.test))
X2.pred[-(nrow(x.test)+1), ] <- x.test
if (semimetric=="deriv") step.2 <- funopare.kNN(y=y, x=x, pred=X2.pred, k=k.opt ,q=m.opt, range.grid=range.grid, nknot=nknot, kind.of.kernel=kind.of.kernel, semimetric=semimetric)
else if (semimetric=="interv") step.2 <- funopare.kNN(y=y, x=x, pred=X2.pred, k=k.opt, interv=c(ext.inf,ext.sup), range.grid=range.grid, nknot=nknot, kind.of.kernel=kind.of.kernel, semimetric=semimetric)
else if (semimetric=="pca") step.2 <- funopare.kNN(y=y, x=x, pred=X2.pred,k=k.opt, q=m.opt, kind.of.kernel=kind.of.kernel, semimetric=semimetric)
Y.estimate.test <- step.2$predicted.values[-(nrow(x.test)+1)]
if(!is.null(y.test))MSE.test.opt <- mean((y.test-Y.estimate.test)^2)
else MSE.test.opt<-NULL
if (semimetric=="interv") {		                                       
	list(m.opt=m.opt, interv.opt=c(ext.inf, ext.sup), k.opt=k.opt, MSE.test.opt=MSE.test.opt, y.estimate.test=Y.estimate.test, 
	CV.app=CV.app, y.test=y.test, parameters.function=parameters.function)                                         
}
else
	list(m.opt=m.opt, k.opt=k.opt, MSE.test.opt=MSE.test.opt, y.estimate.test=Y.estimate.test,
	CV.app=CV.app, y.test=y.test, parameters.function=parameters.function)	
}


