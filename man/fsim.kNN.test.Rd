\name{fsim.kNN.test}
\alias{fsim.kNN.test}

\title{
Functional single-index kNN predictor
}
\description{
Provides predictions when we compute a functional single-index model (FSIM) using the \eqn{k}NN procedure between a scalar response and a functional covariate given a functional index (\eqn{\theta}), a global number of neighbours (\code{k}) and new observations of the functional covariate (\code{x.test}).
}
\usage{
fsim.kNN.test(x, y, x.test, y.test = NULL, theta, order.Bspline = 3,
  nknot.theta = 3, k = 4, kind.of.kernel = "quad", range.grid = NULL, 
  nknot = NULL)
}

\arguments{
  \item{x}{Matrix containing the observations of the functional covariate  that correspond to the training sample collected by row.
}
  \item{y}{Vector containing the scalar responses in the training sample.
}
  \item{x.test}{
Matrix containing the observations of the functional covariate  that correspond to the testing sample collected by row.
}
  \item{y.test}{(optional) Vector/matrix containing the scalar responses in the testing sample.
}
  \item{theta}{
Vector containing the coefficients of \eqn{\theta} in a B-spline basis, so that \code{length(theta)=order.Bspline+nknot.theta}.
}
  \item{order.Bspline}{
Positive integer giving the order of the B-spline basis functions for the B-spline representation of \eqn{\theta}. This is the number of coefficients in each piecewise polynomial segment. The default is 3.
}
  \item{nknot.theta}{
Positive integer indicating the number of uniform interior knots of the B-spline basis for B-spline representation of \eqn{\theta}. The default is 3.
}
  \item{k}{Positive integer indicating the global number of neighbours.}
  \item{kind.of.kernel}{
The type of kernel function used. Only Epanechnikov kernel (\code{"quad"}) is available.
}
  \item{range.grid}{
Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretization). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the size of the discretization size of \code{x} (i.e. \code{ncol(x)}).
}
  \item{nknot}{
Positive integer indicating the number of interior knots for the B-spline representation of the functional covariate. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
}
\details{
The functional single-index model (FSIM) is given by the expression:
\deqn{Y_i=r(\langle\theta_0,X_i\rangle)+\varepsilon_i, \quad i=1,\dots,n,}
where \eqn{Y_i} denotes a scalar response, 
\eqn{X_i} is a functional covariate valued in a separable Hilbert space \eqn{\mathcal{H}} with inner product \eqn{\langle \cdot, \cdot\rangle}, \eqn{\varepsilon} denotes the random error, \eqn{\theta_0 \in \mathcal{H}} is the unknown functional index, \eqn{r(\cdot)} denotes the unknown smooth link function and \eqn{n} is the training sample size.

Given \eqn{\theta \in \mathcal{H}}, \eqn{1<k<n} and a testing sample \{\eqn{X_j,\ j=1,\dots,n_{test}}\}, the predicted responses (see the value \code{y.estimated.test}) can be computed using the kNN procedure by means of 
\deqn{
\widehat{r}_{k,\theta}(X_j)=\sum_{i=1}^nw_{n,k,\theta}(X_j,X_i)Y_i,\quad  j=1,\dots,n_{test}, 
}
with Nadaraya-Watson weights
\deqn{
w_{n,k,\theta}(X_j,X_i)=\frac{K\left(H_{k,X_j,{\theta}}^{-1}d_{\theta}\left(X_i,X_j\right)\right)}{\sum_{i=1}^nK\left(H_{k,X_j,\theta}^{-1}d_{\theta}\left(X_i,X_j\right)\right)},
}

where

\itemize{
\item \eqn{K} is a kernel function (see the argument \code{kind.of.kernel}).

\item for \eqn{x_1,x_2 \in \mathcal{H}, }  \eqn{d_{\theta}(x_1,x_2)=|\langle\theta,x_1-x_2\rangle|} is the projection semi-metric, computed using \code{\link{semimetric.projec}}. \item \eqn{H_{k,x,\theta}=\min\left\{h\in R^+ \text{ such that } \sum_{i=1}^n1_{B_{\theta}(x,h)}(X_i)=k\right\}}, where \eqn{1_{B_{\theta}(x,h)}(\cdot)} is the indicator function of the open ball created with the projection semi-metric with centre \eqn{x\in\mathcal{H}} and radius \eqn{h}.
}
If the argument \code{y.test} is given to the program (i. e. \code{if(!is.null(y.test))}), the function provides the mean squared error of prediction (see the value \code{MSE.test}) calculated as \code{mean((y.test-y.estimated.test)^2)}.

}


\value{
\item{y.estimated.test}{Predicted responses.}
\item{MSE.test}{Mean squared error between predicted and observed responses in the testing sample.}
}

\references{
Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single--index regression. \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
See also \code{\link{fsim.kNN.fit}} and \code{\link{predict.fsim.kNN}}.

Alternative procedure \code{\link{fsim.kernel.test}}.
}
\examples{
\donttest{
data(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra2


train<-1:160
test<-161:215

#FSIM fit. With nknot.theta=2 and range.grid=c(850,1050),
#Theta_n contains 108 thetas.
ptm<-proc.time()
fit<-fsim.kNN.fit(y=y[train],x=X[train,],max.knn=20,nknot.theta=2,nknot=20,
      range.grid=c(850,1050))
proc.time()-ptm
fit

#FSIM prediction
test<-fsim.kNN.test(y=y[train],x=X[train,],x.test=X[test,],y.test=y[test],
        theta=fit$theta.est,k=fit$k.opt,nknot.theta=2,order.Bspline=3,nknot=20,
        range.grid=c(850,1050))

#MSEP
test$MSE.test
}
  
}
