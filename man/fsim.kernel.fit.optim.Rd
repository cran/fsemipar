\name{fsim.kernel.fit.optim}
\alias{fsim.kernel.fit.optim}
\title{
Functional single-index model fit using kernel estimation and iterative LOOCV minimisation
}
\description{
This function fits a functional single-index model (FSIM) between a functional covariate and a
scalar response. 
It employs kernel estimation with Nadaraya-Watson weights and uses B-spline expansions to represent curves and eligible functional indexes. 

The function also utilises the leave-one-out cross-validation (LOOCV) criterion to select the bandwidth (\code{h.opt}) and the coefficients of the functional index in the spline basis (\code{theta.est}). It performs an iterative minimisation of the LOOCV objective function, starting from an initial set of coefficients (\code{gamma}) for the functional index.
}
\usage{
fsim.kernel.fit.optim(x, y, nknot.theta = 3, order.Bspline = 3, gamma = NULL, 
min.q.h = 0.05, max.q.h = 0.5, h.seq = NULL, num.h = 10,
kind.of.kernel = "quad", range.grid = NULL, nknot = NULL, threshold = 0.005)
}
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
  \item{min.q.h}{Minimum quantile order of the distances between curves, which are computed using the projection semi-metric. This value determines the lower endpoint of the range from which the bandwidth is selected. The default is 0.05.
}
  \item{max.q.h}{
Maximum quantile order of the distances between curves, which are computed using the projection semi-metric. This value determines the upper endpoint of the range from which the bandwidth is selected. The default is 0.5.
}
  \item{h.seq}{Vector containing the sequence of bandwidths. The default is a sequence of \code{num.h} equispaced bandwidths in the range constructed using \code{min.q.h} and \code{max.q.h}.
}
  \item{num.h}{Positive integer indicating the number of bandwidths in the grid. The default is 10.
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

The FSIM is fitted using the kernel estimator
\deqn{
\widehat{r}_{h,\hat{\theta}}(x)=\sum_{i=1}^nw_{n,h,\hat{\theta}}(x,X_i)Y_i,  \quad   \forall x\in\mathcal{H},
}
with Nadaraya-Watson weights
\deqn{
w_{n,h,\hat{\theta}}(x,X_i)=\frac{K\left(h^{-1}d_{\hat{\theta}}\left(X_i,x\right)\right)}{\sum_{i=1}^nK\left(h^{-1}d_{\hat{\theta}}\left(X_i,x\right)\right)},
}
 where
 \itemize{
\item the real positive number \eqn{h} is the bandwidth.
 
\item \eqn{K} is a kernel function (see the argument \code{kind.of.kernel}).

\item \eqn{d_{\hat{\theta}}(x_1,x_2)=|\langle\hat{\theta},x_1-x_2\rangle|} is the projection semi-metric, and \eqn{\hat{\theta}} is an estimate of \eqn{\theta_0}. 
}

The procedure requires the estimation of the function-parameter \eqn{\theta_0}. Therefore, we use B-spline expansions to represent curves (dimension \code{nknot+order.Bspline}) and eligible functional indexes (dimension \code{nknot.theta+order.Bspline}).
We obtain the estimated coefficients of \eqn{\theta_0} in the spline basis (\code{theta.est}) and the selected bandwidth (\code{h.opt}) by minimising the LOOCV criterion. This function performs an iterative minimisation procedure, starting from an initial set of coefficients (\code{gamma}) for the functional index. Given a functional index, the optimal bandwidth according to the LOOCV criterion is selected. For a given bandwidth, the minimisation in the functional index is performed using the R function \code{optim}. The procedure is iterated until convergence.  For details, see Ferraty et al. (2013). 
}

\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}.}
\item{theta.est}{Coefficients of \eqn{\hat{\theta}} in the B-spline basis: a vector of \code{length(order.Bspline+nknot.theta)}.}
\item{h.opt}{Selected bandwidth.}
\item{r.squared}{Coefficient of determination.} 
\item{var.res}{Redidual variance.}
\item{df}{Residual degrees of freedom.}
\item{CV.opt}{Minimum value of the LOOCV function, i.e. the value of LOOCV for \code{theta.est} and \code{h.opt}.}
\item{err}{Value of the LOOCV function divided by \code{var(y)} for each interaction.}
\item{H}{Hat matrix.}
\item{h.seq}{Sequence of eligible values for the bandwidth.}
\item{CV.hseq}{CV values for each \code{h}.}
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
See also \code{\link{predict.fsim.kernel}} and \code{\link{plot.fsim.kernel}}.

Alternative procedures \code{\link{fsim.kNN.fit.optim}}, \code{\link{fsim.kernel.fit}} and \code{\link{fsim.kNN.fit}}.
}
\examples{
\donttest{
data(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra2

#FSIM fit.
ptm<-proc.time()
fit<-fsim.kernel.fit.optim(y[1:160],x=X[1:160,],max.q.h=0.35, nknot=20,
range.grid=c(850,1050),nknot.theta=4)
proc.time()-ptm
fit
names(fit)
}
}
