\name{sfpl.kNN.fit}
\alias{sfpl.kNN.fit}
\title{
SFPLM regularised fit using kNN estimation
}
\description{
This function fits a sparse semi-functional partial linear model (SFPLM). It employs a penalised least-squares regularisation procedure, integrated with nonparametric kNN estimation using Nadaraya-Watson weights.

The procedure utilises an objective criterion (\code{criterion}) to select both the bandwidth (\code{h.opt}) and the regularisation parameter (\code{lambda.opt}).
}
\usage{
sfpl.kNN.fit(x, z, y, semimetric = "deriv", q = NULL, knearest = NULL,
min.knn = 2, max.knn = NULL, step = NULL, range.grid = NULL, 
kind.of.kernel = "quad", nknot = NULL, lambda.min = NULL, lambda.min.h = NULL, 
lambda.min.l = NULL, factor.pn = 1, nlambda = 100, lambda.seq = NULL, 
vn = ncol(z), nfolds = 10, seed = 123, criterion = "GCV", penalty = "grSCAD", 
max.iter = 1000)
}

\arguments{
   \item{x}{
Matrix containing the observations of the functional covariate (functional nonparametric component), collected by row.
}
  \item{z}{
Matrix containing the observations of the scalar covariates (linear component), collected by row.
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
  \item{knearest}{
Vector of positive integers containing the sequence in which the  number of nearest neighbours \code{k.opt} is selected. If \code{knearest=NULL}, then \code{knearest <- seq(from =min.knn, to = max.knn, by = step)}.
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
  \item{range.grid}{
Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretisation). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the discretisation size of \code{x} (i.e. \code{ncol(x))}.
}
 \item{kind.of.kernel}{
The type of kernel function used. Currently, only Epanechnikov kernel (\code{"quad"}) is available.
}
  \item{nknot}{
Positive integer indicating the number of interior knots for the B-spline expansion of the functional covariate. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
\item{lambda.min}{
The smallest value for lambda (i.e. the lower endpoint  of the sequence in which \code{lambda.opt} is selected), as fraction of \code{lambda.max}.
The defaults is \code{lambda.min.l} if the sample size is larger than \code{factor.pn} times the number of linear covariates and \code{lambda.min.h} otherwise.
}
  \item{lambda.min.h}{
The lower endpoint of the sequence in which \code{lambda.opt} is selected if the sample size is smaller than \code{factor.pn} times the number of linear covariates. The default is 0.05. 
}
  \item{lambda.min.l}{
The lower endpoint of the sequence in which \code{lambda.opt} is selected if the sample size is larger than \code{factor.pn} times the number of linear covariates. The default is 0.0001.
}
   \item{factor.pn}{
Positive integer used to set \code{lambda.min}. The default value is 1.
  
}
 \item{nlambda}{
 Positive integer indicating the number of values in the sequence from which \code{lambda.opt} is selected. The default is 100.
}
  \item{lambda.seq}{
Sequence of values in which \code{lambda.opt} is selected. If \code{lambda.seq=NULL}, then the programme builds the sequence automatically using \code{lambda.min} and \code{nlambda}.
}
 \item{vn}{
Positive integer or vector of positive integers indicating the number of groups of consecutive variables to be penalised together. The default value is \code{vn=ncol(z)}, resulting in the individual penalization of each scalar covariate.
}
  \item{nfolds}{
Number of cross-validation folds (used when \code{criterion="k-fold-CV"}). Default is 10.
}
  \item{seed}{
You may set the seed for the random number generator to ensure reproducible results (applicable when \code{criterion="k-fold-CV"} is used). The default seed value is 123.
}
   \item{criterion}{
The criterion used to select the tuning and regularisation parameter: \code{k.opt}  and \code{lambda.opt} (also \code{vn.opt} if needed). Options include \code{"GCV"}, \code{"BIC"}, \code{"AIC"}, or \code{"k-fold-CV"}. The default setting is \code{"GCV"}.
}
  \item{penalty}{
The penalty function applied in the penalised least-squares procedure. Currently, only "grLasso" and "grSCAD" are implemented. The default is "grSCAD".
}
  \item{max.iter}{
Maximum number of iterations allowed across the entire path. The default value is 1000.
}
}
\details{
The sparse semi-functional partial linear model (SFPLM) is given by the expression:
\deqn{
Y_i = Z_{i1}\beta_{01} + \dots + Z_{ip_n}\beta_{0p_n} + m(X_i) + \varepsilon_i,\ \ \ i = 1, \dots, n,
}
where \eqn{Y_i} denotes a scalar response, \eqn{Z_{i1}, \dots, Z_{ip_n}} are real random covariates, and \eqn{X_i} is a functional random covariate valued in a semi-metric space \eqn{\mathcal{H}}. In this equation,
\eqn{\mathbf{\beta}_0 = (\beta_{01}, \dots, \beta_{0p_n})^{\top}} and \eqn{m(\cdot)} represent a vector of unknown real parameters and an unknown smooth real-valued function, respectively. Additionally, \eqn{\varepsilon_i} is the random error.

In this function, the SFPLM is fitted using a penalised least-squares approach. The approach involves transforming the SFPLM into a linear model by extracting from \eqn{Y_i} and \eqn{Z_{ij}} (\eqn{j = 1, \ldots, p_n}) the effect of the functional covariate \eqn{X_i} using functional nonparametric regression (for details, see Ferraty and Vieu, 2006). This transformation is achieved using kNN estimation with Nadaraya-Watson weights.

An approximate linear model is then obtained:
\deqn{\widetilde{\mathbf{Y}}\approx\widetilde{\mathbf{Z}}\mathbf{\beta}_0+\mathbf{\varepsilon},}
and the penalised least-squares procedure is applied to this model by minimising
\deqn{
\mathcal{Q}\left(\mathbf{\beta}\right)=\frac{1}{2}\left(\widetilde{\mathbf{Y}}-\widetilde{\mathbf{Z}}\mathbf{\beta}\right)^{\top}\left(\widetilde{\mathbf{Y}}-\widetilde{\mathbf{Z}}\mathbf{\beta}\right)+n\sum_{j=1}^{p_n}\mathcal{P}_{\lambda_{j_n}}\left(|\beta_j|\right), \quad (1)
}
where \eqn{\mathbf{\beta} = (\beta_1, \ldots, \beta_{p_n})^{\top}, \ \mathcal{P}_{\lambda_{j_n}}(\cdot)} is a penalty function (specified in the argument \code{penalty}) and \eqn{\lambda_{j_n} > 0} is a tuning parameter.
To reduce the number of tuning parameters, \eqn{\lambda_j}, to be selected for each sample, we consider \eqn{\lambda_j = \lambda \widehat{\sigma}_{\beta_{0,j,OLS}}}, where \eqn{\beta_{0,j,OLS}} denotes the OLS estimate of \eqn{\beta_{0,j}} and \eqn{\widehat{\sigma}_{\beta_{0,j,OLS}}} is the estimated standard deviation. Both \eqn{\lambda} and \eqn{k} (in the kNN estimation) are selected using the objective criterion specified in the argument \code{criterion}.

Finally, after estimating \eqn{\mathbf{\beta}_0} by minimising (1), we address the estimation of the nonlinear function \eqn{m(\cdot)}.
For this, we again employ the kNN procedure with Nadaraya-Watson weights to smooth the partial residuals \eqn{Y_i - \mathbf{Z}_i^{\top}\widehat{\mathbf{\beta}}}.

For further details on the estimation procedure of the sparse SFPLM, see Aneiros et al. (2015).

\bold{Remark}: It should be noted that if we set \code{lambda.seq} to \eqn{0}, we can obtain the non-penalised estimation of the model, i.e. the OLS estimation. Using \code{lambda.seq} with a value \eqn{\not= 0} is advisable when suspecting the presence of irrelevant variables.
}

\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}}
\item{beta.est}{Estimate of \eqn{\beta_0} when the optimal tuning parameters \code{lambda.opt}, \code{k.opt} and \code{vn.opt} are used.}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{k.opt}{Selected number of nearest neighbours.}
\item{lambda.opt}{Selected value of lambda.}
\item{IC}{Value of the criterion function considered to select both \code{lambda.opt}, \code{h.opt} and \code{vn.opt}.}
\item{vn.opt}{Selected value of \code{vn}.}
\item{...}{}
}
\references{
Aneiros, G., Ferraty, F., Vieu, P. (2015) Variable selection in partial linear regression with functional
covariate. \emph{Statistics}, \bold{49}, 1322--1347, \doi{https://doi.org/10.1080/02331888.2014.998675}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
See also \code{\link{predict.sfpl.kNN}} and  \code{\link{plot.sfpl.kNN}}.

Alternative method \code{\link{sfpl.kernel.fit}}.
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
fit<-sfpl.kNN.fit(y=y[train],x=X[train,], z=z.com[train,],q=2, max.knn=20,
  lambda.min.l=0.01, criterion="BIC",
  range.grid=c(850,1050), nknot=20, max.iter=5000)
proc.time()-ptm

#Results
fit
names(fit)
}
