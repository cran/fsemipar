\name{sfpl.kernel.fit}
\alias{sfpl.kernel.fit}
\title{
Sparse semi-functional partial linear model fit using kernel estimation
}
\description{
This function fits a sparse semi-functional partial linear model between a scalar response, a functional explanatory variable and
a vector of scalar covariates. 
The function uses the penalised least-squares regularization procedure combined with nonparametric kernel estimation with Nadaraya-Watson weights.

The procedure requires an objective criterion (\code{criterion}) to select the bandwidth (\code{h.opt}) and the regularization parameter (\code{lambda.opt}).
}
\usage{
sfpl.kernel.fit(x, z, y, semimetric = "deriv", q = NULL, min.q.h = 0.05,
  max.q.h = 0.5, h.seq = NULL,num.h = 10,range.grid = NULL, 
  kind.of.kernel = "quad", nknot = NULL, lambda.min = NULL, 
  lambda.min.h = NULL, lambda.min.l = NULL, factor.pn = 1, 
  nlambda = 100, lambda.seq = NULL, vn = ncol(z), nfolds = 10, seed = 123,
  criterion = c("GCV", "BIC", "AIC", "k-fold-CV"),
  penalty = c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gBridge", "gLasso",
  "gMCP"),max.iter = 1000)
}
\arguments{
  \item{x}{
Matrix containing the observations of the functional covariate collected by row.
}
  \item{z}{
Matrix containing the observations of the scalar covariates collected by row.
}
  \item{y}{
Vector containing the scalar response.
}
  \item{semimetric}{
Semi-metric function.  Only \code{"deriv"} and \code{"pca"} are implemented.  By default \code{semimetric="deriv"}.
}
  \item{q}{
Order of the derivative (if \code{semimetric="deriv"}) or number of principal components (if \code{semimetric="pca"}). The default values are 0 and 2, respectively.
}
  \item{min.q.h}{
Order of the quantile of the set of distances between curves (computed with the provided \code{semimetric}) which gives the lower end of the sequence in which the bandwidth is selected. The default is 0.05.
}
  \item{max.q.h}{
Order of the quantile of the set of distances between curves (computed with the provided \code{semimetric}) which gives the upper end of the sequence in which the bandwidth is selected. The default is 0.5.
}
  \item{h.seq}{
Vector containing the sequence of bandwidths. The default is a sequence of \code{num.h} equispaced bandwidths in the range constructed using \code{min.q.h} and \code{max.q.h}.
}
  \item{num.h}{
Positive integer indicating the number of bandwiths in the grid. The default is 10.
}
  \item{range.grid}{
Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretization). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the size of the discretization size of \code{x} (i.e. \code{ncol(x))}.
}
  \item{kind.of.kernel}{
The type of kernel function used. Only Epanechnikov kernel (\code{"quad"}) is available.
}
  \item{nknot}{
Positive integer indicating the number of interior knots for the B-spline representation of the functional covariate. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
  \item{lambda.min}{
The smallest value for lambda (i. e., the smallest value  of the sequence in which \code{lambda.opt} is selected), as fraction of \code{lambda.max}.
The defaults is \code{lambda.min.l} if the number of observations is larger than \code{factor.pn} times the number of covariates and \code{lambda.min.h} otherwise.
}
  \item{lambda.min.h}{
The smallest value of the sequence in which \code{lambda.opt} is selected if the number of observations is smaller than \code{factor.pn} times the number of scalar covariates. The default is 0.05. 
}
  \item{lambda.min.l}{
  The smallest value of the sequence in which \code{lambda.opt} is selected if the number of observations is larger than \code{factor.pn} times the number of scalar covariates. The default is 0.0001.
}
  \item{factor.pn}{
  Positive integer used to set \code{lambda.min}. The default is 1.
}
  \item{nlambda}{
  Positive integer indicating the number of values of the sequence in which \code{lambda.opt} is selected. The default is 100.
}
  \item{lambda.seq}{
Sequence of values in which \code{lambda.opt} is selected. If \code{lambda.seq=NULL}, then the programme builds the sequence automatically using \code{lambda.min} and \code{nlambda}.
}
  \item{vn}{
Positive integer or vector of positive integers indicating the number of groups of consecutive variables to be penalised together. The default value is \code{vn=ncol(z)}, which leads to the individual penalisation of each scalar covariate.
}
  \item{nfolds}{
Positive integer indicating the number of cross-validation folds (used if \code{criterion="k-fold-CV"}). The default is 10.
}
  \item{seed}{
You may set the seed of the random number generator to obtain reproducible results (used if \code{criterion="k-fold-CV"}). The default is 123.
}
  \item{criterion}{
The criterion by which to select the regularization parameter \code{lambda.opt} and the bandwidth \code{h.opt}. One of \code{"GCV", "BIC", "AIC"} or \code{"k-fold-CV"}. The default is \code{"GCV"}.
}
  \item{penalty}{
The penalty function to be applied in the penalized least squares procedure. Only "grLasso" and "grSCAD" are implemented.
}

\item{max.iter}{
 Maximum number of iterations (total across entire path). The default is 1000.
}

}
\details{
The sparse semi-functional partial linear model (SSFPLM) is given by the expression:
\deqn{
Y_i=Z_{i1}\beta_{01}+\dots+Z_{ip_n}\beta_{0p_n}+m(X_i)+\varepsilon_i\ \ \ i=1,\dots,n,
}
where \eqn{Y_i} denotes a scalar response, \eqn{Z_{i1},\dots,Z_{ip_n}} are real random covariates and \eqn{X_i} is a functional random covariate valued in some semi-metric space \eqn{\mathcal{H}}. In this equation,
	\eqn{\mathbf{\beta}_0=(\beta_{01},\dots,\beta_{0p_n})^{\top}} and \eqn{m(\cdot)} are a vector of unknown real parameters and an unknown smooth real-valued function, respectively. In addition, \eqn{\varepsilon_i} is the random error.


In this function the SSFPLM is fitted using the penalised least-squares approach. The first idea is to transform the SSFPLM into a linear model by extracting from \eqn{Y_i} and \eqn{Z_{ij}} (\eqn{j=1,\ldots,p_n}) the effect of the functional covariate \eqn{X_i} using functional nonparametric regression (see, for details, Ferraty and Vieu, 2006). This is made using  kernel estimation with Nadaraya-Watson weights.

Then, an approximate linear model is obtained:
\deqn{\widetilde{\mathbf{Y}}\approx\widetilde{\mathbf{Z}}\mathbf{\beta}_0+\mathbf{\varepsilon},}
and the penalised least-squares procedure is applied to this model by minimising 
\deqn{
\mathcal{Q}\left(\mathbf{\beta}\right)=\frac{1}{2}\left(\widetilde{\mathbf{Y}}-\widetilde{\mathbf{Z}}\mathbf{\beta}\right)^{\top}\left(\widetilde{\mathbf{Y}}-\widetilde{\mathbf{Z}}\mathbf{\beta}\right)+n\sum_{j=1}^{p_n}\mathcal{P}_{\lambda_{j_n}}\left(|\beta_j|\right), \quad (1)
}

where \eqn{\mathbf{\beta}=(\beta_1,\ldots,\beta_{p_n})^{\top}, \ \mathcal{P}_{\lambda_{j_n}}\left(\cdot\right)} is a penalty function (specified in the argument \code{penalty}) and \eqn{\lambda_{j_n} > 0} is a tuning parameter.
To reduce  the quantity of tuning parameters, \eqn{\lambda_j}, to be selected for each sample, we consider \eqn{\lambda_j = \lambda \widehat{\sigma}_{\beta_{0,j,OLS}}}, where \eqn{\beta_{0,j,OLS}} denotes the OLS estimate of \eqn{\beta_{0,j}} and \eqn{\widehat{\sigma}_{\beta_{0,j,OLS}}} is the estimated standard deviation. Both \eqn{\lambda} and \eqn{h} (in the kernel estimation) are selected using the objetive criterion specified in the argument \code{criterion}.


Finally, after estimating \eqn{\mathbf{\beta}_0} by minimising (1), we deal with the estimation of the nonlinear function \eqn{m(\cdot)}.
For that, we employ again the  kernel procedure with Nadaraya-Watson weights to smooth the partial residuals \eqn{Y_i-\mathbf{Z}_i^{\top}\widehat{\mathbf{\beta}}}.

For further details on the estimation procedure of the SSFPLM, see Aneiros et al. (2015).

\bold{Remark}: We should note that if we set \code{lambda.seq}\eqn{=0}, we can obtain the non-penalised estimation of the model, i.e. the OLS estimation. It is convenient to use \code{lambda.seq}\eqn{\not=0} when one suspects there are irrelevant variables.

}
\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}}
\item{beta.est}{Estimate of \eqn{\beta_0} when the optimal tuning parameters \code{lambda.opt}, \code{h.opt} and \code{vn.opt} are used.}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{h.opt}{Selected bandwidth.}
\item{lambda.opt}{Selected value of lambda.}
\item{IC}{Value of the criterion function considered to select \code{lambda.opt}, \code{h.opt} and \code{vn.opt}.}
\item{h.min.opt.max.mopt}{\code{h.opt=h.min.opt.max.mopt[2]} (used by \code{beta.est}) was seeked between \code{h.min.opt.max.mopt[1]} and \code{h.min.opt.max.mopt[3]}.}
\item{vn.opt}{Selected value of \code{vn}.}
\item{...}{}
}
\references{
Aneiros, G., Ferraty, F., Vieu, P. (2015) Variable selection in partial linear regression with functional
covariate. \emph{Statistics}, \bold{49}, 1322--1347, \doi{https://doi.org/10.1080/02331888.2014.998675}.

Ferraty, F. and Vieu, P. (2006) \emph{Nonparametric Functional Data Analysis}. Springer Series in Statistics, New York.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}

\seealso{
See also \code{\link{predict.sfpl.kernel}} and  \code{\link{plot.sfpl.kernel}}.

Alternative method \code{\link{sfpl.kNN.fit}}.
}
\examples{
data("Tecator")
y<-Tecator$fat
X<-Tecator$absor.spectra
z1<-Tecator$protein       
z2<-Tecator$moisture

#Quadratic, cubic and interaction effects of the scalar covariates.
z.com<-cbind(z1,z2,z1^2,z2^2,z1^3,z2^3,z1*z2)
train<-1:160


#SFPLM fit. 
ptm=proc.time()
fit<-sfpl.kernel.fit(x=X[train,], z=z.com[train,], y=y[train],q=2, 
      max.q.h=0.35,lambda.min.h=0.02,lambda.min.l=0.01,
      factor.pn=2, max.iter=5000, criterion="BIC", penalty="grSCAD",nknot=20)
proc.time()-ptm

#Results
fit
names(fit)
}

