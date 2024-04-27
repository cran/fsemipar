\name{predict.sfpl}
\alias{predict.sfpl.kernel}
\alias{predict.sfpl.kNN}
\title{
Predictions for SFPLM
}
\description{
\code{predict} method for the semi-functional partial linear  model (SFPLM) fitted using \code{sfpl.kernel.fit} or \code{sfpl.kNN.fit}. 
}

\usage{
\method{predict}{sfpl.kernel}(object, newdata.x = NULL, newdata.z = NULL,
  y.test = NULL, option = NULL, ...)
\method{predict}{sfpl.kNN}(object, newdata.x = NULL, newdata.z = NULL, 
  y.test = NULL, option = NULL, ...)
}

\arguments{
  \item{object}{
Output of the functions mentioned in the \code{Description} (i.e. an object of the class \code{sfpl.kernel} or \code{sfpl.kNN}.
}
  \item{newdata.x}{
Matrix containing new observations of the functional covariate collected by row.
}
  \item{newdata.z}{
Matrix containing the new observations of the scalar covariate collected by row.
}
  \item{y.test}{
(optional) A vector containing the new observations of the response.
}
  \item{option}{
Allows the selection among the choices 1 and 2 for \code{sfpl.kernel} objects, and 1, 2 and 3 for \code{sfpl.kNN} objects. The default setting is 1. See the section \code{Details}.
}
  \item{\dots}{Further arguments passed to or from other methods.
}
}
\details{
 The following options are provided to obtain the predictions of the response for \code{newdata.x} and \code{newdata.z}: 
\itemize{
\item If \code{option=1}, we maintain all the estimations (\code{k.opt} or \code{h.opt} and \code{beta.est}) to predict the functional nonparametric component of the model. 
\item If \code{option=2}, we maintain \code{beta.est}, while the tuning parameter (\eqn{h} or \eqn{k}) is selected again to predict the functional nonparametric component of the model. This selection is performed using the leave-one-out cross-validation (LOOCV) criterion in the associated functional nonparametric model, obtaining a global selection for \eqn{h} or \eqn{k}.
}
In the case of \code{sfpl.kNN} objects if \code{option=3}, we maintain \code{beta.est}, while the tuning parameter \eqn{k} is seleted again to predict the functional nonparametric component of the model. This selection is performed using the LOOCV criterion in the associated functional nonparametric model, performing a local selection for \eqn{k}.
}
\value{
The function returns the predicted values of the response (\code{y}) for \code{newdata.x} and \code{newdata.z}. If \code{!is.null(y.test)}, it also provides the mean squared error of prediction (\code{MSEP}) computed as \code{mean((y-y.test)^2)}.
If \code{is.null(newdata.x)} or \code{is.null(newdata.z)}, then the function returns the fitted values.
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}

\seealso{
\code{\link{sfpl.kernel.fit}} and \code{\link{sfpl.kNN.fit}}
}
\examples{
\donttest{
data("Tecator")
y<-Tecator$fat
X<-Tecator$absor.spectra
z1<-Tecator$protein       
z2<-Tecator$moisture

#Quadratic, cubic and interaction effects of the scalar covariates.
z.com<-cbind(z1,z2,z1^2,z2^2,z1^3,z2^3,z1*z2)
train<-1:160
test<-161:215

 
#Fit
fit.kernel<-sfpl.kernel.fit(x=X[train,], z=z.com[train,], y=y[train],q=2,
  max.q.h=0.35,lambda.min.l=0.01, factor.pn=2,  
  criterion="BIC", range.grid=c(850,1050), nknot=20, max.iter=5000)
fit.kNN<-sfpl.kNN.fit(y=y[train],x=X[train,], z=z.com[train,],q=2, 
  max.knn=20,lambda.min.l=0.01, factor.pn=2, 
  criterion="BIC",range.grid=c(850,1050), nknot=20, max.iter=5000)

#Predictions
predict(fit.kernel,newdata.x=X[test,],newdata.z=z.com[test,],y.test=y[test],
  option=2)
predict(fit.kNN,newdata.x=X[test,],newdata.z=z.com[test,],y.test=y[test],
  option=2)
}
}
