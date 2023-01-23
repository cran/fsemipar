\name{fsim.kNN.fit}
\alias{fsim.kNN.fit}
\title{Functional single-index model fit using kNN estimation}
\description{
This function fits a functional single-index model (FSIM) between a functional explanatory variable and
scalar response. 
The function uses \eqn{k}-nearest neighbours (kNN) estimation with Nadaraya-Watson weights, a  B-spline representation to estimate the functional index \eqn{\theta_0} and the cross-validation (CV) criterion to select the number of neighbours (\code{k.opt}) and the coefficients of the functional index in the spline basis (\code{theta.est}).
}
\usage{
fsim.kNN.fit(x, y, seed.coeff = c(-1, 0, 1), order.Bspline = 3, 
  nknot.theta = 3, t0 = NULL, min.knn = 2, max.knn = NULL, knearest = NULL,
  step = NULL, kind.of.kernel = "quad", range.grid = NULL, nknot = NULL)
}

\arguments{
  \item{x}{Matrix containing the observations of the functional covariate collected by row.}
  
  \item{y}{Vector containing the scalar response.}
  
  \item{seed.coeff}{Vector of initial values used to  build the set \eqn{\Theta_n} (see section \bold{Details}). The coefficients for the B-spline representation of each eligible functional index \eqn{\theta \in \Theta_n} are obtained from \code{seed.coeff}.  The default is \code{c(-1,0,1)}.}
  
  \item{order.Bspline}{Positive integer giving the order of the B-spline basis functions. This is the number of coefficients in each piecewise polynomial segment. The default is 3.}
  
  \item{nknot.theta}{Positive integer indicating the number of uniform interior knots of the B-spline basis for the B-spline representation of \eqn{\theta_0}. The default is 3.}
  
  \item{t0}{Value in the domain of the functional indexes at which we evaluate them to  build the set \eqn{\Theta_n}. We assume \eqn{\theta_0(t_0)>0} for some arbitrary \eqn{t_0} in the domain to ensure model identifiability. If \code{t0=NULL}, then \code{mean(range.grid)} is considered.}
  
  \item{min.knn}{Positive integer indicating the smallest value of the sequence in which the  number of nearest neighbours \code{k.opt} is selected (thus, this number must be smaller than the sample size). The default is 2.}
  
  \item{max.knn}{Positive integer indicating the largest value of the sequence in which the number of nearest neighbours \code{k.opt} is selected (thus, this number must be larger than \code{min.kNN} and smaller than the sample size, \code{n}). The default is \code{max.knn <- n\%/\%2}.
}
  \item{knearest}{Vector of positive integers containing the sequence in which the  number of nearest neighbours \code{k.opt} is selected. If \code{knearest=NULL}, then \code{knearest <- seq(from =min.knn, to = max.knn, by = step)}.}
  
  \item{step}{Positive integer used to build the sequence of k-nearest neighbours in the following way: \code{min.knn, min.knn + step, min.knn + 2*step, min.knn + 3*step,...} The default is \code{step<-ceiling(n/100)}.}
  
  \item{kind.of.kernel}{The type of kernel function used. Only Epanechnikov kernel (\code{"quad"}) is available.
              
}
  \item{range.grid}{Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretization). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the size of the discretization size of \code{x} (i.e. \code{ncol(x))}.}
  
  \item{nknot}{Positive integer indicating the number of interior knots for the B-spline representation of the functional covariate. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
}

\details{

The functional single-index model (FSIM) is given by the expression:
\deqn{Y_i=r(\langle\theta_0,X_i\rangle)+\varepsilon_i, \quad i=1,\dots,n,}
where \eqn{Y_i} denotes a scalar response, 
\eqn{X_i} is a functional covariate valued in a separable Hilbert space \eqn{\mathcal{H}} with inner product \eqn{\langle \cdot, \cdot\rangle}, \eqn{\varepsilon} denotes the random error, \eqn{\theta_0 \in \mathcal{H}} is the unknown functional index and \eqn{r(\cdot)} denotes the unknown smooth link function.

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
\item \eqn{k\in \mathbb{Z}^+} is a smoothing factor (the number of nearest neighbours).
 
\item \eqn{K} is a kernel function (see the argument \code{kind.of.kernel}).

\item \eqn{d_{\hat{\theta}}(x_1,x_2)=|\langle\hat{\theta},x_1-x_2\rangle|} is the projection semi-metric, computed using \code{\link{semimetric.projec}} and \eqn{\hat{\theta}} is an estimate of \eqn{\theta_0}. 

\item \eqn{H_{k,x,\hat{\theta}}=\min\{h\in \mathbb{R}^+ \text{ such that } \sum_{i=1}^n1_{B_{\hat{\theta}}(x,h)}(X_i)=k\}}, where \eqn{1_{B_{\hat{\theta}}(x,h)}(\cdot)} is the indicator function of the open ball created with the projection semi-metric with centre \eqn{x\in\mathcal{H}} and radius \eqn{h}.
}

The procedure requires the estimation of the function-parameter \eqn{\theta_0}. Therefore, we use B-spline representation to build a set  \eqn{\Theta_n} of eligible functional indexes. The dimension of the B-spline basis is \code{order.Bspline}+\code{nknot.theta} and the set of eligible coefficients is obtained by calibrating (to ensure the identifiability of the model) the set of initial coefficients given in \code{seed.coeff}. The larger this set, the higher the size of \eqn{\Theta_n}. Since our approach requires intensive computation, we need a trade-off between the size of \eqn{\Theta_n} and the performance of the estimator. For that, Ait-Saidi et al. (2008) suggested considering \code{order.Bspline=3} and \code{seed.coeff=c(-1,0,1)}. For details on the construction of \eqn{\Theta_n} see Novo et al. (2019).

We obtain the estimated coefficients of \eqn{\theta_0} in the spline basis (\code{theta.est}) and the selected number of neighbours (\code{k.opt}) by minimising the CV criterion.

}


\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}}
\item{theta.est}{Coefficients of \eqn{\hat{\theta}} in the B-spline basis: a vector of \code{length(order.Bspline+nknot.theta)}.}
\item{k.opt}{Selected number of nearest neighbours.}
\item{r.squared}{Coefficient of determination.} 
\item{var.res}{Redidual variance.}
\item{df}{Residual degrees of freedom.}
\item{yhat.cv}{Predicted values for the scalar response using leave-one-out samples.}
\item{CV.opt}{Minimum value of the CV function, i.e. the value of CV for \code{theta.est} and \code{k.opt}.}
\item{CV.values}{Vector containing CV values for each functional index in \eqn{\Theta_n} and the value of \eqn{k} that minimises the CV for such index (i.e. \code{CV.values[j]} contains the value of the CV function corresponding to \code{theta.seq.norm[j,]} and the best value of the \code{k.seq} for this functional index according to the CV criterion).}
\item{H}{Hat matrix.}
\item{m.opt}{Index of \eqn{\hat{\theta}} in the set \eqn{\Theta_n}.}
\item{theta.seq.norm}{The vector \code{theta.seq.norm[j,]} contains the coefficientes in the B-spline basis of the jth functional index in \eqn{\Theta_n}.}
\item{k.seq}{Sequence of eligible values for \eqn{k}.}
\item{...}{}
}


\references{
Ait-Saidi, A., Ferraty, F., Kassa, R., and Vieu, P. (2008) Cross-validated estimations in the single-functional index model, \emph{Statistics}, \bold{42(6)}, 475--494, \doi{https://doi.org/10.1080/02331880801980377}.

Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single--index regression, \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}




\seealso{
See also \code{\link{fsim.kNN.test}}, \code{\link{predict.fsim.kNN}}, \code{\link{plot.fsim.kNN}}.

Alternative procedure \code{\link{fsim.kernel.fit}}.
}
\examples{
\donttest{

data(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra2

#FSIM fit. With nknot.theta=2 and range.grid=c(850,1050), 
#Theta_n contains 108 thetas.
ptm<-proc.time()
fit<-fsim.kNN.fit(y=y[1:160],x=X[1:160,],max.knn=20,nknot.theta=2,nknot=20,
range.grid=c(850,1050))
proc.time()-ptm
fit
names(fit)
}

}



