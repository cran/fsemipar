\name{fsim.kernel.fit}
\alias{fsim.kernel.fit}

\title{
Functional single-index model fit using kernel estimation
}
\description{
This function fits a functional single-index model (FSIM) between a functional explanatory variable and
scalar response. 
The function uses kernel estimation with Nadaraya-Watson weights, a  B-spline representation to estimate the functional index \eqn{\theta_0} and the cross-validation (CV) criterion to select the bandwidth (\code{h.opt}) and the coefficients of the functional index in the spline basis (\code{theta.est}).
}
\usage{
fsim.kernel.fit(x, y, seed.coeff = c(-1, 0, 1), nknot.theta = 3, 
  order.Bspline = 3, t0 = NULL, min.q.h = 0.05, max.q.h = 0.5,
  h.seq = NULL, num.h = 10, kind.of.kernel = "quad", range.grid = NULL,
  nknot = NULL)
}

\arguments{
  \item{x}{Matrix containing the observations of the functional covariate (i.e. curves) collected by row.}
  
  \item{y}{Vector containing the scalar response.}
  
  \item{seed.coeff}{
Vector of initial values used to  build the set \eqn{\Theta_n} (see section \code{Details}). The coefficients for the B-spline representation of each eligible functional index \eqn{\theta \in \Theta_n} are obtained from \code{seed.coeff}.  The default is \code{c(-1,0,1)}.
}
  \item{nknot.theta}{Positive integer indicating the number of uniform interior knots of the B-spline basis for the B-spline representation of \eqn{\theta_0}. The default is 3.
}
  \item{order.Bspline}{
Positive integer giving the order of the B-spline basis functions. This is the number of coefficients in each piecewise polynomial segment. The default is 3.
}
  \item{t0}{
Value in the domain of the functional indexes at which we evaluate them to build the set \eqn{\Theta_n}. We assume \eqn{\theta_0(t_0)>0} for some arbitrary \eqn{t_0} in the domain to ensure model identifiability. If \code{t0=NULL}, then \code{mean(range.grid)} is considered.
}
  \item{min.q.h}{Order of the quantile of the set of distances between curves (computed with  the projection semi-metric) which gives the lower end of the sequence in which the bandwidth is selected. The default is 0.05.
}
  \item{max.q.h}{
Order of the quantile of the set of distances between curves (computed with the projection semi-metric) which gives the upper end of the sequence in which the bandwidth is selected. The default is 0.5.
}
  \item{h.seq}{Vector containing the sequence of bandwidths. The default is a sequence of \code{num.h} equispaced bandwidths in the range constructed using \code{min.q.h} and \code{max.q.h}.
}
  \item{num.h}{Positive integer indicating the number of bandwiths in the grid. The default is 10.
}
  \item{kind.of.kernel}{
The type of kernel function used. Only Epanechnikov kernel (\code{"quad"}) is available.
}
  \item{range.grid}{
Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretization). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the size of the discretization size of \code{x} (i.e. \code{ncol(x))}.
}
  \item{nknot}{
Positive integer indicating the number of interior knots for the B-spline representation of the functional covariate. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
}
\details{
The functional single-index model (FSIM) is given by the expression:
\deqn{Y_i=r(\langle\theta_0,X_i\rangle)+\varepsilon_i, \quad i=1,\dots,n,}
where \eqn{Y_i} denotes a scalar response, 
\eqn{X_i} is a functional covariate valued in a separable Hilbert space \eqn{\mathcal{H}} with inner product \eqn{\langle \cdot, \cdot\rangle}, \eqn{\varepsilon} denotes the random error, \eqn{\theta_0 \in \mathcal{H}} is the unknown functional index and \eqn{r(\cdot)} denotes the unknown smooth link function.

The FSIM is fitted using the kNN estimator
\deqn{
\widehat{r}_{h,\hat{\theta}}(x)=\sum_{i=1}^nw_{n,h,\hat{\theta}}(x,X_i)Y_i,  \quad   \forall x\in\mathcal{H},
}
with Nadaraya-Watson weights
\deqn{
w_{n,h,\hat{\theta}}(x,X_i)=\frac{K\left(h^{-1}d_{\hat{\theta}}\left(X_i,x\right)\right)}{\sum_{i=1}^nK\left(h^{-1}d_{\hat{\theta}}\left(X_i,x\right)\right)},
}
 where
 \itemize{
\item \eqn{h\in \mathbb{R}^+} is the bandwidth.
 
\item \eqn{K} is a kernel function (see the argument \code{kind.of.kernel}).

\item \eqn{d_{\hat{\theta}}(x_1,x_2)=|\langle\hat{\theta},x_1-x_2\rangle|} is the projection semi-metric, computed using \code{\link{semimetric.projec}} and \eqn{\hat{\theta}} is an estimate of \eqn{\theta_0}. 
}

The procedure requires the estimation of the function-parameter \eqn{\theta_0}. Therefore, we use B-spline representation to build a set  \eqn{\Theta_n} of eligible functional indexes. The dimension of the B-spline basis is \code{order.Bspline}+\code{nknot.theta} and the set of eligible coefficients is obtained by calibrating (to ensure the identifiability of the model) the set of initial coefficients given in \code{seed.coeff}. The larger this set, the higher the size of \eqn{\Theta_n}. Since our approach requires intensive computation, we need a trade-off between the size of \eqn{\Theta_n} and the performance of the estimator. For that, Ait-Saidi et al. (2008) suggested considering \code{order.Bspline=3} and \code{seed.coeff=c(-1,0,1)}. For details on the construction of \eqn{\Theta_n} see Novo et al. (2019).

We obtain the estimated coefficients of \eqn{\theta_0} in the spline basis (\code{theta.est}) and the selected bandwidth (\code{h.opt}) by minimising the CV criterion.

}
\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}}
\item{theta.est}{Coefficients of \eqn{\hat{\theta}} in the B-spline basis: a vector of \code{length(order.Bspline+nknot.theta)}.}
\item{h.opt}{Selected bandwidth.}
\item{r.squared}{Coefficient of determination.} 
\item{var.res}{Redidual variance.}
\item{df}{Residual degrees of freedom.}
\item{yhat.cv}{Predicted values for the scalar response using leave-one-out samples.}
\item{CV.opt}{Minimum value of the CV function, i.e. the value of CV for \code{theta.est} and \code{h.opt}.}
\item{CV.values}{Vector containing CV values for each functional index in \eqn{\Theta_n} and the value of \eqn{h} that minimises the CV for such index (i.e. \code{CV.values[j]} contains the value of the CV function corresponding to \code{theta.seq.norm[j,]} and the best value of the \code{h.seq} for this functional index according to the CV criterion).}
\item{H}{Hat matrix.}
\item{m.opt}{Index of \eqn{\hat{\theta}} in the set \eqn{\Theta_n}.}
\item{theta.seq.norm}{The vector \code{theta.seq.norm[j,]} contains the coefficientes in the B-spline basis of the jth functional index in \eqn{\Theta_n}.}
\item{h.seq}{Sequence of eligible values for \eqn{h}.}
\item{...}{}
}
\references{
Ait-Saidi, A., Ferraty, F., Kassa, R., and Vieu, P. (2008) Cross-validated estimations in the single-functional index model. \emph{Statistics}, \bold{42(6)}, 475--494, \doi{https://doi.org/10.1080/02331880801980377}.

Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single--index regression. \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
See also \code{\link{fsim.kernel.test}}, \code{\link{predict.fsim.kernel}}, \code{\link{plot.fsim.kernel}}.

Alternative procedure \code{\link{fsim.kNN.fit}}.
}
\examples{
\donttest{
data(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra2


#FSIM fit. With nknot.theta=2 and range.grid=c(850,1050), 
#Theta_n contains 108 thetas.
ptm<-proc.time()
fit<-fsim.kernel.fit(y[1:160],x=X[1:160,],max.q.h=0.35, nknot=20,
range.grid=c(850,1050),nknot.theta=2)
proc.time()-ptm
fit
names(fit)
}
}


