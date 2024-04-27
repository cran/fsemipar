\name{predict.mfplm.PVS}
\alias{predict.PVS.kernel}
\alias{predict.PVS.kNN}
\title{
Prediction for MFPLM
}
\description{
\code{predict} method for the multi-functional partial linear model (MFPLM) fitted using \code{PVS.kernel.fit} or \code{PVS.kNN.fit}.
}

\usage{
\method{predict}{PVS.kernel}(object, newdata.x = NULL, newdata.z = NULL,
  y.test = NULL, option = NULL, ...)
\method{predict}{PVS.kNN}(object, newdata.x = NULL, newdata.z = NULL, 
  y.test = NULL, option = NULL, knearest.n = object$knearest, 
  min.knn.n = object$min.knn, max.knn.n = object$max.knn.n, 
  step.n = object$step, ...)

}
\arguments{
   \item{object}{
Output of the functions mentioned in the \code{Description} (i.e. an object of the class \code{PVS.kernel} or \code{PVS.kNN}).
}
  \item{newdata.x}{
A matrix containing new observations of the functional covariate in the functional nonparametric component, collected by row.
}
  \item{newdata.z}{
Matrix containing the new observations of the scalar covariates derived from the discretisation  of a curve,  collected by row. 
}
  \item{y.test}{
(optional) A vector containing the new observations of the response.
}
  \item{option}{
Allows the selection among the choices 1, 2 and 3 for \code{PVS.kernel} objects, and  1, 2, 3, and 4 for \code{PVS.kNN} objects. The default setting is 1. See the section \code{Details}.
}

 \item{\dots}{
Further arguments.
}
   \item{knearest.n}{
Only used for objects \code{PVS.kNN} if \code{option=2}, \code{option=3} or \code{option=4}: sequence in which the  number of nearest neighbours \code{k.opt} is selected. The default is \code{object$knearest}.
}
  \item{min.knn.n}{
Only used for objects \code{PVS.kNN} if \code{option=2}, \code{option=3} or \code{option=4}: minumum value of the sequence in which the  number of nearest neighbours \code{k.opt} is selected (thus, this number must be smaller than the sample size). The default is \code{object$min.knn}.
}
  \item{max.knn.n}{
Only used for objects \code{PVS.kNN} if \code{option=2}, \code{option=3} or \code{option=4}: maximum value of the sequence in which the number of nearest neighbours \code{k.opt} is selected (thus, this number must be larger than \code{min.kNN} and smaller than the sample size). The default is \code{object$max.knn}.
}
  \item{step.n}{
Only used for objects \code{PVS.kNN} if \code{option=2}, \code{option=3} or \code{option=4}: positive integer used to build the sequence of k-nearest neighbours in the following way: \code{min.knn, min.knn + step.n, min.knn + 2*step.n, min.knn + 3*step.n,...}. The default is  \code{object$step}.
}
}
\details{
To obtain the predictions of the response for \code{newdata.x} and \code{newdata.z}, the following options are provided:
\itemize{
\item If \code{option=1}, we maintain all the estimates (\code{k.opt} or \code{h.opt} and \code{beta.est}) to predict the functional nonparametric component of the model. As we use the estimates of the second step of the algorithm, only the \code{train.2} is used as training sample to predict.
Then, it should be noted that \code{k.opt} or \code{h.opt} may not be suitable to predict the functional nonparametric component of the model.
\item If \code{option=2}, we maintain \code{beta.est}, while the tuning parameter (\eqn{h} or \eqn{k}) is selected again to predict the functional nonparametric component of the model. This selection is performed using the leave-one-out cross-validation (LOOCV) criterion in the associated functional nonparametric model and the complete training sample (i.e. \code{train=c(train.1,train.2)}), obtaining a global selection for \eqn{h} or \eqn{k}. As we use the entire training sample (not just a subsample of it), the sample size is modified and, as a consequence,  the parameters \code{knearest}, \code{min.knn}, \code{max.knn}, and \code{step} given to the function \code{IASSMR.kNN.fit} may need to be provided again to compute predictions. For that, we add the arguments \code{knearest.n}, \code{min.knn.n}, \code{max.knn.n} and \code{step.mn}. 
\item  If \code{option=3}, we maintain only the indexes of the relevant variables selected by the IASSMR. We estimate again the linear coefficients  using \code{sfpl.kernel.fit} or \code{sfpl.kNN.fit}, respectively, without penalisation (setting \code{lambda.seq=0}) and using the entire training sample (\code{train=c(train.1,train.2)}). The method provides two predictions (and MSEPs):
\itemize{
\item a) The prediction associated with \code{option=1} for \code{sfpl.kernel} or \code{sfpl.kNN} class.
\item b) The prediction associated with \code{option=2} for \code{sfpl.kernel} or \code{sfpl.kNN} class.
}
(see the documentation of the functions \code{predict.sfpl.kernel} and \code{predict.sfpl.kNN})
\item If \code{option=4} (an option only available for the class \code{PVS.kNN}) we maintain \code{beta.est}, while the tuning parameter \eqn{k} is selected again to predict the functional nonparametric component of the model. This selection is performed using LOOCV criterion in the functional nonparametric model associated and the complete training sample (i.e. \code{train=c(train.1,train.2)}), obtaining a local selection for \eqn{k}.
}
}
\value{
The function returns the predicted values of the response (\code{y}) for \code{newdata.x} and \code{newdata.z}. If \code{!is.null(y.test)}, it also provides the mean squared error of prediction (\code{MSEP}) computed as \code{mean((y-y.test)^2)}.
If \code{option=3}, two sets of predictions (and two MSEPs) are provided, corresponding to the items a) and b) mentioned in the section \code{Details.}
If \code{is.null(newdata.x)} or \code{is.null(newdata.z)}, then the function returns the fitted values.
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}



\seealso{
\code{\link{PVS.kernel.fit}}, \code{\link{sfpl.kernel.fit}} and \code{\link{predict.sfpl.kernel}} or \code{\link{PVS.kNN.fit}},
 \code{\link{sfpl.kNN.fit}} and \code{\link{predict.sfpl.kNN}}.
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
fit.kernel<- PVS.kernel.fit(x=x.sug[train,],z=z.sug[train,], 
              y=y.sug[train],train.1=1:108,train.2=109:216,
              lambda.min.h=0.03,lambda.min.l=0.03,
              max.q.h=0.35, nknot=20,criterion="BIC",
              max.iter=5000)
fit.kNN<- PVS.kNN.fit(x=x.sug[train,],z=z.sug[train,], y=y.sug[train],
            train.1=1:108,train.2=109:216,lambda.min.h=0.07, 
            lambda.min.l=0.07, nknot=20,criterion="BIC",
            max.iter=5000)

#Preditions
predict(fit.kernel,newdata.x=x.sug[test,],newdata.z=z.sug[test,],y.test=y.sug[test],option=2)
predict(fit.kNN,newdata.x=x.sug[test,],newdata.z=z.sug[test,],y.test=y.sug[test],option=2)
}
}
