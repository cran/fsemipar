\name{predict.lm}
\alias{predict.PVS}
\alias{predict.lm.pels}
\title{
Prediction from linear model estimates
}
\description{
\code{predict} method for: 
\itemize{
\item Linear model (LM) fitted using \code{lm.pels.fit}. 
\item Linear model with covariates coming from the discretization of a curve fitted using \code{PVS.fit}.
}
}
\usage{
\method{predict}{lm.pels}(object, newdata = NULL, y.test = NULL, ...)
\method{predict}{PVS}(object, newdata = NULL, y.test = NULL, ...)
}
\arguments{
  \item{object}{
Output of the \code{lm.pels.fit} or \code{PVS.fit} functions (i.e. an object of the class \code{lm.pels} or \code{PVS})
}
  \item{newdata}{
Matrix containing the new observations of the scalar covariates (LM) or of the scalar covariates coming from the discretisation  of a curve, collected by row.
}
  \item{y.test}{
(optional) A vector containing the new observations of the response.
}
  \item{\dots}{
}
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
\code{\link{lm.pels.fit}} and \code{\link{PVS.fit}}.
}
\examples{

data("Tecator")
y<-Tecator$fat
z1<-Tecator$protein       
z2<-Tecator$moisture

#Quadratic, cubic and interaction effects of the scalar covariates.
z.com<-cbind(z1,z2,z1^2,z2^2,z1^3,z2^3,z1*z2)
train<-1:160
test<-161:215

#LM fit. 
fit<-lm.pels.fit(z=z.com[train,], y=y[train],lambda.min.h=0.02,lambda.min.l=0.01,
      factor.pn=2, max.iter=5000, criterion="BIC", penalty="grSCAD")

#Predictions
predict(fit,newdata=z.com[test,],y.test=y[test])


data(Sugar)

y<-Sugar$ash
z<-Sugar$wave.240

#Outliers
index.y.25 <- y > 25
index.atip <- index.y.25
(1:268)[index.atip]


#Dataset to model
z.sug<- z[!index.atip,]
y.sug <- y[!index.atip]

train<-1:216
test<-217:266

#Fit
fit.pvs<-PVS.fit(z=z.sug[train,], y=y.sug[train],train.1=1:108,train.2=109:216,
          lambda.min.h=0.2,criterion="BIC", penalty="grSCAD", max.iter=5000)


#Predictions
predict(fit.pvs,newdata=z.sug[test,],y.test=y.sug[test])


}

