\name{predict.fsim}

\alias{predict.fsim.kernel}
\alias{predict.fsim.kNN}

\title{
Prediction for FSIM
}
\description{
\code{predict} method for the functional single-index model (FSIM) fitted using \code{fsim.kernel.fit}, \code{fsim.kernel.fit.optim}, \code{fsim.kNN.fit} and \code{fsim.kNN.fit.optim}.
}
\usage{
\method{predict}{fsim.kernel}(object, newdata = NULL, y.test = NULL, ...)
\method{predict}{fsim.kNN}(object, newdata = NULL, y.test = NULL, ...)
}

\arguments{
  \item{object}{
Output of the \code{fsim.kernel.fit}, \code{fsim.kernel.fit.optim}, \code{fsim.kNN.fit} or \code{fsim.kNN.fit.optim} functions (i.e. an object of the class \code{fsim.kernel} or \code{fsim.kNN}).
}
  \item{newdata}{
A matrix containing new observations of the functional covariate collected by row.
}
  \item{y.test}{
(optional) A vector containing the new observations of the response.
}
  \item{\dots}{Further arguments passed to or from other methods.

}
}
\details{
The prediction is computed using the functions \code{fsim.kernel.test} and \code{fsim.kernel.fit}, respectively.
}
\value{
The function returns the predicted values of the response (\code{y}) for \code{newdata}. If \code{!is.null(y.test)}, it also provides the mean squared error of prediction (\code{MSEP}) computed as \code{mean((y-y.test)^2)}.
If \code{is.null(newdata)} the function returns the fitted values.
}


\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}

\seealso{
\code{fsim.kernel.fit} and  \code{fsim.kernel.test}  or \code{fsim.kNN.fit} and \code{fsim.kNN.test}.
}
\examples{
\donttest{
data(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra2

train<-1:160
test<-161:215

#FSIM fit. 
fit.kernel<-fsim.kernel.fit(y[train],x=X[train,],max.q.h=0.35, nknot=20,
range.grid=c(850,1050),nknot.theta=4)
fit.kNN<-fsim.kNN.fit(y=y[train],x=X[train,],max.knn=20,nknot=20,
nknot.theta=4, range.grid=c(850,1050))

test<-161:215

pred.kernel<-predict(fit.kernel,newdata=X[test,],y.test=y[test])
pred.kernel$MSEP
pred.kNN<-predict(fit.kNN,newdata=X[test,],y.test=y[test])
pred.kNN$MSEP
}
}

