\name{predict.sfplsim.FASSMR}
\alias{predict.FASSMR.kernel}
\alias{predict.FASSMR.kNN}
\alias{predict.sfplsim.kernel}
\alias{predict.sfplsim.kNN}

\title{
Prediction from functional semiparametric  model estimates
}
\description{
\code{predict} method for: 
\itemize{
\item Semi-functional partial linear single-index model (SFPLSIM) fitted using \code{sfplsim.kernel.fit} or \code{sfplsim.kNN.fit}. 
\item Multi-functional partial linear single-index regression (MFPLSIM) fitted using \code{FASSMR.kernel.fit} or \code{FASSMR.kNN.fit}.
}
}
\usage{
\method{predict}{sfplsim.kernel}(object, newdata.x = NULL, newdata.z = NULL,
  y.test = NULL, option = NULL, ...)
\method{predict}{sfplsim.kNN}(object, newdata.x = NULL, newdata.z = NULL, 
  y.test = NULL, option = NULL, ...)
\method{predict}{FASSMR.kernel}(object, newdata.x = NULL, newdata.z = NULL,
  y.test = NULL, option = NULL, ...)
\method{predict}{FASSMR.kNN}(object, newdata.x = NULL, newdata.z = NULL,
  y.test = NULL, option = NULL, ...)
}

\arguments{
  \item{object}{
Output of the functions mentioned in the \code{Description} (i.e. an object of the class \code{sfplsim.kernel}, \code{sfplsim.kNN}, \code{FASSMR.kernel} or \code{FASSMR.kNN}).
}
  \item{newdata.x}{
A matrix containing new observations of the functional covariate in the functional-single index component collected by row.
}
  \item{newdata.z}{
Matrix containing the new observations of the scalar covariates (SFPLSIM) or of the scalar covariates coming from the discretisation  of a curve (MFPLSIM), collected by row. 
}
  \item{y.test}{
(optional) A vector containing the new observations of the response.
}
  \item{option}{
Allows the choice between 1 and 2. The default is 1. See the section \code{Details}.
}
  \item{\dots}{Further arguments passed to or from other methods.

}
}
\details{
To obtain the predictions of the response for \code{newdata.x} and \code{newdata.z}, two options are provided:
\itemize{
\item If \code{option=1}, we mantain all the estimations (\code{k.opt} or \code{h.opt}, \code{theta.est} and \code{beta.est}) to predict the functional single-index component of the model. 
\item If \code{option=2}, we mantain \code{theta.est} and \code{beta.est}, while the tuning parameter (\eqn{h} or \eqn{k}) is selected again to predict the functional single-index component of the model. This selection is performed using the cross-validation criterion in the functional single-index model associated.
}
}
\value{
The function returns the predicted values of the response (\code{y}) for \code{newdata.x} and \code{newdata.z}. If \code{!is.null(y.test)}, it also provides the mean squared error of prediction (\code{MSEP}) computed as \code{mean((y-y.test)^2)}.
If \code{is.null(newdata.x)} or \code{is.null(newdata.z)},  the function returns the fitted values.
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}

\seealso{
\code{\link{sfplsim.kernel.fit}}, \code{\link{sfplsim.kNN.fit}}, \code{\link{FASSMR.kernel.fit}}  or \code{\link{FASSMR.kNN.fit}}.
}
\examples{
\donttest{
data("Tecator")
y<-Tecator$fat
X<-Tecator$absor.spectra2
z1<-Tecator$protein       
z2<-Tecator$moisture

#Quadratic, cubic and interaction effects of the scalar covariates.
z.com<-cbind(z1,z2,z1^2,z2^2,z1^3,z2^3,z1*z2)
train<-1:160
test<-161:215

#SFPLSIM fit. Convergence errors for some theta are obtained.
s.fit.kernel<-sfplsim.kernel.fit(x=X[train,], z=z.com[train,], y=y[train],
            max.q.h=0.35,lambda.min.h=0.02,lambda.min.l=0.01, 
            factor.pn=2, max.iter=5000, nknot.theta=4,criterion="BIC", 
            penalty="grSCAD",nknot=20)
s.fit.kNN<-sfplsim.kNN.fit(y=y[train],x=X[train,], z=z.com[train,], 
        max.knn=20,lambda.min.h=0.02,lambda.min.l=0.01, factor.pn=2, 
        nknot.theta=4,criterion="BIC",range.grid=c(850,1050), penalty="grSCAD",
        nknot=20, max.iter=5000)


predict(s.fit.kernel,newdata.x=X[test,],newdata.z=z.com[test,],
  y.test=y[test],option=2)
predict(s.fit.kNN,newdata.x=X[test,],newdata.z=z.com[test,],
  y.test=y[test],option=2)
}
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

m.fit.kernel <- FASSMR.kernel.fit(x=x.sug[train,],z=z.sug[train,], 
                  y=y.sug[train],  nknot.theta=2,lambda.min.h=0.03, 
                  lambda.min.l=0.03, max.q.h=0.35,num.h = 10, 
                  nknot=20,criterion="BIC", penalty="grSCAD",max.iter=5000)


m.fit.kNN<- FASSMR.kNN.fit(x=x.sug[train,],z=z.sug[train,], y=y.sug[train], 
            nknot.theta=2,lambda.min.h=0.03, lambda.min.l=0.03, 
            max.knn=20,nknot=20,criterion="BIC", penalty="grSCAD",max.iter=5000)


predict(m.fit.kernel,newdata.x=x.sug[test,],newdata.z=z.sug[test,],
  y.test=y.sug[test],option=2)
predict(m.fit.kNN,newdata.x=x.sug[test,],newdata.z=z.sug[test,],
  y.test=y.sug[test],option=2)
}


}

