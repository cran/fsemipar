\name{fsim.kNN.fit.optim}
\alias{fsim.kNN.fit.optim}

\title{
Functional single-index model fit using kNN estimation and iterative LOOCV minimisation
}
\description{
This function fits a functional single-index model (FSIM) between a functional covariate and a
scalar response. 
It employs kNN estimation with Nadaraya-Watson weights and uses B-spline expansions to represent curves and eligible functional indexes. 

The function also utilises the leave-one-out cross-validation (LOOCV) criterion to select the bandwidth (\code{h.opt}) and the coefficients of the functional index in the spline basis (\code{theta.est}). It performs an iterative minimisation of the LOOCV objective function, starting from an initial set of coefficients (\code{gamma}) for the functional index.
}
\usage{
fsim.kNN.fit.optim(x, y, order.Bspline = 3, nknot.theta = 3, gamma = NULL, 
knearest = NULL, min.knn = 2, max.knn = NULL,  step = NULL, 
kind.of.kernel = "quad", range.grid = NULL, nknot = NULL, threshold = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{Matrix containing the observations of the functional covariate (i.e. curves) collected by row.}
  
  \item{y}{Vector containing the scalar response.}
  
  \item{order.Bspline}{
Positive integer giving the order of the B-spline basis functions. This is the number of coefficients in each piecewise polynomial segment. The default is 3}

  \item{nknot.theta}{Positive integer indicating the number of regularly spaced interior knots in the B-spline expansion of \eqn{\theta_0}. The default is 3. 
}
  \item{gamma}{
Vector indicating the initial coefficients for the functional index used in the iterative procedure. By default, it is a vector of ones. The size of the vector is determined by the sum \code{nknot.theta+order.Bspline}.
}

  \item{knearest}{
Vector of positive integers that defines the sequence within which the optimal number of nearest neighbours \code{k.opt} is selected. If \code{knearest=NULL}, then \code{knearest <- seq(from =min.knn, to = max.knn, by = step)}.
}
  \item{min.knn}{
A positive integer that represents the minimum value in the sequence for selecting the number of nearest neighbours \code{k.opt}. This value should be less than the sample size. The default is 2.
}
  \item{max.knn}{
A positive integer that represents the maximum value in the sequence for selecting number of nearest neighbours \code{k.opt}. This value should be less than the sample size. The default is \code{max.knn <- n\%/\%5}.
}
   \item{step}{
A positive integer used to construct the sequence of k-nearest neighbours as follows: \code{min.knn, min.knn + step, min.knn + 2*step, min.knn + 3*step,...}. The default value for \code{step} is \code{step<-ceiling(n/100)}.
}
   \item{kind.of.kernel}{
The type of kernel function used. Currently, only Epanechnikov kernel (\code{"quad"}) is available.
}
 \item{range.grid}{
Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretisation). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the discretisation size of \code{x} (i.e. \code{ncol(x))}.
}
  \item{nknot}{
Positive integer indicating the number of regularly spaced interior knots for the B-spline expansion of the functional covariate. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
  \item{threshold}{
The convergence threshold for the LOOCV function (scaled by the variance of the response). The default is \code{5e-3}.
}
}
\details{
The functional single-index model (FSIM) is given by the expression:
\deqn{Y_i=r(\langle\theta_0,X_i\rangle)+\varepsilon_i, \quad i=1,\dots,n,}
where \eqn{Y_i} denotes a scalar response, 
\eqn{X_i} is a functional covariate valued in a separable Hilbert space \eqn{\mathcal{H}} with an inner product \eqn{\langle \cdot, \cdot\rangle}. The term \eqn{\varepsilon} denotes the random error, \eqn{\theta_0 \in \mathcal{H}} is the unknown functional index and \eqn{r(\cdot)} denotes the unknown smooth link function.

The FSIM is fitted using the kNN estimator
\deqn{
\widehat{r}_{k,\hat{\theta}}(x)=\sum_{i=1}^nw_{n,k,\hat{\theta}}(x,X_i)Y_i,  \quad   \forall x\in\mathcal{H},
}
with Nadaraya-Watson weights
\deqn{
w_{n,k,\hat{\theta}}(x,X_i)=\frac{K\left(H_{k,x,\hat{\theta}}^{-1}d_{\hat{\theta}}\left(X_i,x\right)\right)}{\sum_{i=1}^nK\left(H_{k,x,\hat{\theta}}^{-1}d_{\hat{\theta}}\left(X_i,x\right)\right)},
}
 where
 \itemize{
\item the positive integer \eqn{k} is a smoothing factor, representing the number of nearest neighbours.
 
\item \eqn{K} is a kernel function (see the argument \code{kind.of.kernel}).

\item \eqn{d_{\hat{\theta}}(x_1,x_2)=|\langle\hat{\theta},x_1-x_2\rangle|} is the projection semi-metric and \eqn{\hat{\theta}} is an estimate of \eqn{\theta_0}. 

\item \eqn{H_{k,x,\hat{\theta}}=\min\{h\in R^+ \text{ such that } \sum_{i=1}^n1_{B_{\hat{\theta}}(x,h)}(X_i)=k\}}, where \eqn{1_{B_{\hat{\theta}}(x,h)}(\cdot)} is the indicator function of the open ball defined by the projection semi-metric, with centre \eqn{x\in\mathcal{H}} and radius \eqn{h}.
}

The procedure requires the estimation of the function-parameter \eqn{\theta_0}. Therefore, we use B-spline expansions to represent curves (dimension \code{nknot+order.Bspline}) and eligible functional indexes (dimension \code{nknot.theta+order.Bspline}).
We obtain the estimated coefficients of \eqn{\theta_0} in the spline basis (\code{theta.est}) and the selected number of neighbours (\code{k.opt}) by minimising the LOOCV criterion. This function performs an iterative minimisation procedure, starting from an initial set of coefficients (\code{gamma}) for the functional index. Given a functional index, the optimal number of neighbours according to the LOOCV criterion is selected. For a given number of neighbours, the minimisation in the functional index is performed using the R function \code{optim}. The procedure is iterated until convergence.  For details, see Ferraty et al. (2013). 

}
\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}.}
\item{theta.est}{Coefficients of \eqn{\hat{\theta}} in the B-spline basis: a vector of \code{length(order.Bspline+nknot.theta)}.}
\item{k.opt}{Selected number of neighbours.}
\item{r.squared}{Coefficient of determination.} 
\item{var.res}{Redidual variance.}
\item{df}{Residual degrees of freedom.}
\item{CV.opt}{Minimum value of the LOOCV function, i.e. the value of LOOCV for \code{theta.est} and \code{k.opt}.}
\item{err}{Value of the LOOCV function divided by \code{var(y)} for each interaction.}
\item{H}{Hat matrix.}
\item{k.seq}{Sequence of eligible values for \eqn{k}.}
\item{CV.hseq}{CV values for each \eqn{k}.}
\item{...}{}
}
\references{
Ferraty, F., Goia, A., Salinelli, E., and Vieu, P. (2013) Functional projection pursuit regression. \emph{Test}, \bold{22}, 293--320, \doi{https://doi.org/10.1007/s11749-012-0306-2}.

Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single--index regression. \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}
\seealso{
See also \code{\link{predict.fsim.kNN}} and \code{\link{plot.fsim.kNN}}.

Alternative procedures \code{\link{fsim.kernel.fit.optim}}, \code{\link{fsim.kernel.fit}} and \code{\link{fsim.kNN.fit}}.
}
\examples{
\donttest{
data(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra2

#FSIM fit.
ptm<-proc.time()
fit<-fsim.kNN.fit.optim(y=y[1:160],x=X[1:160,],max.knn=20,nknot.theta=4,nknot=20,
range.grid=c(850,1050))
proc.time()-ptm
fit
names(fit)
}

}