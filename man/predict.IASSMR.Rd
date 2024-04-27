\name{predict.IASSMR}
\alias{predict.IASSMR.kernel}
\alias{predict.IASSMR.kNN}

\title{
Prediction for MFPLSIM
}
\description{
\code{predict} method for the multi-functional partial linear single-index model (MFPLSIM) fitted using \code{IASSMR.kernel.fit} or \code{IASSMR.kNN.fit}.
}
\usage{

\method{predict}{IASSMR.kernel}(object, newdata.x = NULL, newdata.z = NULL,
  y.test = NULL, option = NULL, ...)
\method{predict}{IASSMR.kNN}(object, newdata.x = NULL, newdata.z = NULL,
  y.test = NULL, option = NULL, knearest.n = object$knearest, 
  min.knn.n = object$min.knn, max.knn.n = object$max.knn.n, 
  step.n = object$step, ...)

}

\arguments{
  \item{object}{
Output of the functions mentioned in the \code{Description} (i.e. an object of the class \code{IASSMR.kernel} or \code{IASSMR.kNN}).
}
  \item{newdata.x}{
A matrix containing new observations of the functional covariate in the functional single-index component, collected by row.
}
  \item{newdata.z}{
Matrix containing the new observations of the scalar covariates derived from the discretisation  of a curve,  collected by row. 
}
  \item{y.test}{
(optional) A vector containing the new observations of the response.
}
  \item{option}{
Allows the choice between 1, 2 and 3. The default is 1. See the section \code{Details}.
}

 \item{\dots}{
Further arguments.
}
   \item{knearest.n}{
Only used for objects \code{IASSMR.kNN} if \code{option=2} or \code{option=3}: vector of positive integers containing the sequence in which the  number of nearest neighbours \code{k.opt} is selected. The default is \code{object$knearest}.
}
  \item{min.knn.n}{
Only used for objects \code{IASSMR.kNN} if \code{option=2} or \code{option=3}: minumum value of the sequence in which the  number of neighbours \code{k.opt} is selected (thus, this number must be smaller than the sample size). The default is \code{object$min.knn}.
}
  \item{max.knn.n}{
Only used for objects \code{IASSMR.kNN} if \code{option=2} or \code{option=3}: maximum value of the sequence in which the number of neighbours \code{k.opt} is selected (thus, this number must be larger than \code{min.kNN} and smaller than the sample size). The default is \code{object$max.knn}.
}
  \item{step.n}{
Only used for objects \code{IASSMR.kNN} if \code{option=2} or \code{option=3}: positive integer used to build the sequence of k-nearest neighbours as follows: \code{min.knn, min.knn + step.n, min.knn + 2*step.n, min.knn + 3*step.n,...}. The default is  \code{object$step}.
}
 
}
\details{
Three options are provided to obtain the predictions of the response for \code{newdata.x} and \code{newdata.z}:
\itemize{
\item If \code{option=1}, we maintain all the estimates (\code{k.opt} or \code{h.opt}, \code{theta.est} and \code{beta.est}) to predict the functional single-index component of the model. As we use the estimates of the second step of the algorithm, only the \code{train.2} is used as training sample to predict.
Then, it should be noted that \code{k.opt} or \code{h.opt} may not be suitable to predict the functional single-index component of the model.
\item If \code{option=2}, we maintain \code{theta.est} and \code{beta.est}, while the tuning parameter (\eqn{h} or \eqn{k}) is selected again to predict the functional single-index component of the model. This selection is performed using the leave-one-out cross-validation criterion in the functional single-index model associated and the complete training sample (i.e. \code{train=c(train.1,train.2)}). As we use the entire training sample (not just a subsample of it), the sample size is modified and, as a consequence,  the parameters \code{knearest}, \code{min.knn}, \code{max.knn}, \code{step} given to the function \code{IASSMR.kNN.fit} may need to be provided again to compute predictions. For that, we add the arguments \code{knearest.n}, \code{min.knn.n}, \code{max.knn.n} and \code{step.n}.
\item  If \code{option=3}, we maintain only the indexes of the relevant variables selected by the IASSMR. We estimate again the linear coefficients and the functional index  by means of \code{sfplsim.kernel.fit} or \code{sfplsim.kNN.fit}, respectively, without penalisation (setting \code{lambda.seq=0}) and using the whole training sample (\code{train=c(train.1,train.2)}). The method provides two predictions (and MSEPs):
\itemize{
\item a) The prediction associated with \code{option=1} for \code{sfplsim.kernel} or \code{sfplsim.kNN} class.
\item b) The prediction associated with \code{option=2} for \code{sfplsim.kernel} or \code{sfplsim.kNN} class.
}
(see the documentation of the functions \code{predict.sfplsim.kernel} and \code{predict.sfplsim.kNN})

}
}
\value{
The function returns the predicted values of the response (\code{y}) for \code{newdata.x} and \code{newdata.z}. If \code{!is.null(y.test)}, it also provides the mean squared error of prediction (\code{MSEP}) computed as \code{mean((y-y.test)^2)}.
If \code{option=3}, two sets of predictions (and two MSEPs) are provided, corresponding to the items a) and b) mentioned in the section \code{Details.}
If \code{is.null(newdata.x)} or \code{is.null(newdata.z)},  the function returns the fitted values.

}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}

\seealso{
\code{\link{sfplsim.kernel.fit}}, \code{\link{sfplsim.kNN.fit}}, \code{\link{IASSMR.kernel.fit}}, \code{\link{IASSMR.kNN.fit}}.
}
\examples{
\donttest{
data(Sugar)

y<-Sugar$ash
x<-Sugar$wave.290
z<-Sugar$wave.240

#Outliers
index.y.25 <- y > 25
index.atip <- index.y.25
(1:268)[index.atip]

#Dataset to model
x.sug <- x[!index.atip,]
z.sug<- z[!index.atip,]
y.sug <- y[!index.atip]

train<-1:216
test<-217:266

#Fit
fit.kernel<-IASSMR.kernel.fit(x=x.sug[train,],z=z.sug[train,], y=y.sug[train],
            train.1=1:108,train.2=109:216,nknot.theta=2,lambda.min.h=0.03,
            lambda.min.l=0.03,  max.q.h=0.35,  nknot=20,criterion="BIC",
            max.iter=5000)

fit.kNN<- IASSMR.kNN.fit(x=x.sug[train,],z=z.sug[train,], y=y.sug[train],
          train.1=1:108,train.2=109:216,nknot.theta=2,lambda.min.h=0.07,
          lambda.min.l=0.07, max.knn=20, nknot=20,criterion="BIC",
          max.iter=5000)

#Predictions
predict(fit.kernel,newdata.x=x.sug[test,],newdata.z=z.sug[test,],y.test=y.sug[test],option=2)
predict(fit.kNN,newdata.x=x.sug[test,],newdata.z=z.sug[test,],y.test=y.sug[test],option=2)
}
}

