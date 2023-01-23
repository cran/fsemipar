fnp.kernel.fit.test <- function(y,x,x.test,y.test=NULL, 
kind.of.semimetric="semimetric.pca", start.order.deriv.o.pca=NULL, end.order.deriv.o.pca=NULL, min.leng.interv=NULL, max.leng.interv=NULL,
min.quantile.h=0.05, max.quantile.h=0.5, h.seq=NULL, num.h=NULL,  
kind.of.kernel="quad",range.grid=NULL, nknot=NULL)
{
# The construction of this code is based on that by F. Ferraty, which is available on his website 
# https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html.
if (!is.matrix(x))  x <- t(as.matrix(x))
if (!is.matrix(x.test)) x.test <- t(as.matrix(x.test))
ext.inf <- ext.sup <- NULL
aux.cv <- fnp.kernel.fit(x=x,y=y, kind.of.semimetric=kind.of.semimetric, kind.of.kernel=kind.of.kernel, min.quantile.h=min.quantile.h, max.quantile.h=max.quantile.h, h.seq=h.seq, num.h=num.h, 
			start.order.deriv.o.pca=start.order.deriv.o.pca, end.order.deriv.o.pca=end.order.deriv.o.pca, min.leng.interv=min.leng.interv, max.leng.interv=max.leng.interv, range.grid=range.grid, nknot=nknot)
m.opt <- aux.cv$m.opt
h.opt <- aux.cv$h.opt
h.seq.opt <- aux.cv$h.seq.opt
Y.estimate.app <- aux.cv$Y.estimate.app
CV.app <- aux.cv$CV.app
if (kind.of.semimetric=="semimetric.interv") {	
	ext.inf <- aux.cv$interv.opt[1]
	ext.sup <- aux.cv$interv.opt[2]
}
aux.pred <- fnp.kernel.test(x=x,y=y,x.test=x.test,y.test=y.test, kind.of.semimetric=kind.of.semimetric, q=m.opt, h=h.opt, ext.inf=ext.inf, ext.sup=ext.sup, kind.of.kernel=kind.of.kernel, range.grid=range.grid, nknot=nknot)
MSEP <- aux.pred$MSE.test
Y.estimate.test <- aux.pred$y.estimate.test
if (kind.of.semimetric=="semimetric.interv") 		                                          
	list(m.opt=m.opt, h.opt=h.opt, h.seq.opt=h.seq.opt, interv.opt=c(ext.inf, ext.sup), CV=CV.app, MSEP=MSEP, y.estimate.test = Y.estimate.test)                                         
else
	list(m.opt=m.opt, h.opt=h.opt, h.seq.opt=h.seq.opt, CV.app=CV.app, MSEP=MSEP, y.estimate.test = Y.estimate.test)	
}





