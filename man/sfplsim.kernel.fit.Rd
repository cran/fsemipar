\name{sfplsim.kernel.fit}
\alias{sfplsim.kernel.fit}

\title{
SFPLSIM regularised fit using kernel estimation
}
\description{
This function fits a sparse semi-functional partial linear single-index (SFPLSIM). It employs a penalised least-squares regularisation procedure, integrated with nonparametric kernel estimation using Nadaraya-Watson weights.

The function uses B-spline expansions to represent curves and eligible functional indexes.  It also utilises an objective criterion (\code{criterion}) to select both the bandwidth (\code{h.opt}) and the regularisation parameter (\code{lambda.opt}).
}
\usage{
sfplsim.kernel.fit(x, z, y, seed.coeff = c(-1, 0, 1), order.Bspline = 3, 
nknot.theta = 3, min.q.h = 0.05, max.q.h = 0.5, h.seq = NULL, num.h = 10, 
range.grid = NULL, kind.of.kernel = "quad", nknot = NULL, lambda.min = NULL,
lambda.min.h = NULL, lambda.min.l = NULL, factor.pn = 1, nlambda = 100, 
lambda.seq = NULL, vn = ncol(z), nfolds = 10, seed = 123, criterion = "GCV",
penalty = "grSCAD", max.iter = 1000, n.core = NULL)
}

\arguments{
  \item{x}{
Matrix containing the observations of the functional covariate (functional single-index component), collected by row.
}
  \item{z}{
Matrix containing the observations of the scalar covariates (linear component), collected by row.
}
  \item{y}{
Vector containing the scalar response.
}
\item{seed.coeff}{
Vector of initial values used to  build the set \eqn{\Theta_n} (see section \code{Details}). The coefficients for the B-spline representation of each eligible functional index \eqn{\theta \in \Theta_n} are obtained from \code{seed.coeff}.  The default is \code{c(-1,0,1)}.
}
  \item{order.Bspline}{
Positive integer giving the order of the B-spline basis functions. This is the number of coefficients in each piecewise polynomial segment. The default is 3.
}
 \item{nknot.theta}{
Positive integer indicating the number of regularly spaced interior knots in the B-spline expansion of \eqn{\theta_0}. The default is 3.
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
The smallest value for lambda (i. e., the lower endpoint  of the sequence in which \code{lambda.opt} is selected), as fraction of \code{lambda.max}.
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
The criterion used to select the tuning and regularisation parameter: \code{h.opt} and \code{lambda.opt}  (also \code{vn.opt} if needed). Options include \code{"GCV"}, \code{"BIC"}, \code{"AIC"}, or \code{"k-fold-CV"}. The default setting is \code{"GCV"}.
}
  \item{penalty}{
The penalty function applied in the penalised least-squares procedure. Currently, only "grLasso" and "grSCAD" are implemented. The default is "grSCAD".
}
  \item{max.iter}{
Maximum number of iterations allowed across the entire path. The default value is 1000.
}
  \item{n.core}{
Number of CPU cores designated for parallel execution. The default is \code{n.core<-availableCores(omit=1)}.
}
}
\details{
The sparse semi-functional partial linear single-index model (SFPLSIM) is given by the expression:
\deqn{
Y_i=Z_{i1}\beta_{01}+\dots+Z_{ip_n}\beta_{0p_n}+r(\left<\theta_0,X_i\right>)+\varepsilon_i\ \ \ i=1,\dots,n,
}
where \eqn{Y_i} denotes a scalar response, \eqn{Z_{i1},\dots,Z_{ip_n}} are real random covariates and \eqn{X_i} is a functional random covariate valued in a separable Hilbert space \eqn{\mathcal{H}} with inner product \eqn{\left\langle \cdot, \cdot \right\rangle}. In this equation,
	\eqn{\mathbf{\beta}_0=(\beta_{01},\dots,\beta_{0p_n})^{\top}}, \eqn{\theta_0\in\mathcal{H}} and \eqn{r(\cdot)} are a vector of unknown real parameters, an unknown functional direction and an unknown smooth real-valued function, respectively. In addition, \eqn{\varepsilon_i} is the random error.


The sparse SFPLSIM is fitted using the penalised least-squares approach. The first step is to transform the SSFPLSIM into a linear model by extracting from \eqn{Y_i} and \eqn{Z_{ij}} (\eqn{j=1,\ldots,p_n}) the effect of the functional covariate \eqn{X_i} using functional single-index regression.  This transformation is achieved using nonparametric kernel estimation (see, for details, the documentation of the function \code{fsim.kernel.fit}).

An approximate linear model is then obtained:
\deqn{\widetilde{\mathbf{Y}}_{\theta_0}\approx\widetilde{\mathbf{Z}}_{\theta_0}\mathbf{\beta}_0+\mathbf{\varepsilon},}
and the penalised least-squares procedure is applied to this model by minimising over the pair \eqn{(\mathbf{\beta},\theta)}
\deqn{
\mathcal{Q}\left(\mathbf{\beta},\theta\right)=\frac{1}{2}\left(\widetilde{\mathbf{Y}}_{\theta}-\widetilde{\mathbf{Z}}_{\theta}\mathbf{\beta}\right)^{\top}\left(\widetilde{\mathbf{Y}}_{\theta}-\widetilde{\mathbf{Z}}_{\theta}\mathbf{\beta}\right)+n\sum_{j=1}^{p_n}\mathcal{P}_{\lambda_{j_n}}\left(|\beta_j|\right), \quad (1)
}
where \eqn{\mathbf{\beta}=(\beta_1,\ldots,\beta_{p_n})^{\top}, \ \mathcal{P}_{\lambda_{j_n}}\left(\cdot\right)} is a penalty function (specified in the argument \code{penalty}) and \eqn{\lambda_{j_n} > 0} is a tuning parameter.
To reduce  the quantity of tuning parameters, \eqn{\lambda_j}, to be selected for each sample, we consider \eqn{\lambda_j = \lambda \widehat{\sigma}_{\beta_{0,j,OLS}}}, where \eqn{\beta_{0,j,OLS}} denotes the OLS estimate of \eqn{\beta_{0,j}} and \eqn{\widehat{\sigma}_{\beta_{0,j,OLS}}} is the estimated standard deviation. Both \eqn{\lambda} and \eqn{h} (in the kernel estimation) are selected using the objetive criterion specified in the argument \code{criterion}.
 

In addition, the function uses a B-spline representation to construct a set  \eqn{\Theta_n} of eligible functional indexes \eqn{\theta}. The dimension of the B-spline basis is \code{order.Bspline}+\code{nknot.theta} and the set of eligible coefficients is obtained by calibrating (to ensure the identifiability of the model) the set of initial coefficients given in \code{seed.coeff}. The larger this set, the greater the size of \eqn{\Theta_n}. ue to the intensive computation required by our approach, a balance between the size of \eqn{\Theta_n} and the performance of the estimator is necessary. For that, Ait-Saidi et al. (2008) suggested considering \code{order.Bspline=3} and \code{seed.coeff=c(-1,0,1)}. For details on the construction of \eqn{\Theta_n} see Novo et al. (2019).


Finally, after estimating \eqn{\mathbf{\beta}_0} and \eqn{\theta_0} by minimising (1), we proceed to estimate the nonlinear function \eqn{r_{\theta_0}(\cdot)\equiv r\left(\left<\theta_0,\cdot\right>\right)}.
For this purporse, we again apply the kernel procedure with Nadaraya-Watson weights to smooth the partial residuals \eqn{Y_i-\mathbf{Z}_i^{\top}\widehat{\mathbf{\beta}}}.

For further details on the estimation procedure of the SSFPLSIM, see Novo et al. (2021).

\bold{Remark}: It should be noted that if we set \code{lambda.seq} to \eqn{0}, we can obtain the non-penalised estimation of the model, i.e. the OLS estimation. Using \code{lambda.seq} with a value \eqn{\not= 0} is advisable when suspecting the presence of irrelevant variables.
}
\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}.}
\item{beta.est}{Estimate of \eqn{\beta_0} when the optimal tuning parameters \code{lambda.opt}, \code{h.opt} and \code{vn.opt} are used.}
\item{theta.est}{Coefficients of \eqn{\hat{\theta}} in the B-spline basis (when the optimal tuning parameters \code{lambda.opt}, \code{h.opt} and \code{vn.opt} are used): a vector of \code{length(order.Bspline+nknot.theta)}.}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{h.opt}{Selected bandwidth.}
\item{lambda.opt}{Selected value of the penalisation parameter \eqn{\lambda}.}
\item{IC}{Value of the criterion function considered to select \code{lambda.opt}, \code{h.opt} and \code{vn.opt}.}
\item{Q.opt}{Minimum value of the penalized criterion used to estimate \eqn{\beta_0} and \eqn{\theta_0}. That is, the value obtained using \code{theta.est} and \code{beta.est}.}
\item{Q}{Vector of dimension equal to the cardinal of \eqn{\Theta_n}, containing the values of the penalized criterion for each functional index in \eqn{\Theta_n}.}
\item{m.opt}{Index of \eqn{\hat{\theta}} in the set \eqn{\Theta_n}.}
\item{lambda.min.opt.max.mopt}{A grid of values in [\code{lambda.min.opt.max.mopt[1], lambda.min.opt.max.mopt[3]}] is considered to seek for the \code{lambda.opt} (\code{lambda.opt=lambda.min.opt.max.mopt[2]}).
}
\item{lambda.min.opt.max.m}{A grid of values in [\code{lambda.min.opt.max.m[m,1], lambda.min.opt.max.m[m,3]}] is considered to seek for the optimal \eqn{\lambda} (\code{lambda.min.opt.max.m[m,2]})
used by the optimal \eqn{\beta} for each \eqn{\theta} in \eqn{\Theta_n}.}
\item{h.min.opt.max.mopt}{\code{h.opt=h.min.opt.max.mopt[2]} (used by \code{theta.est} and \code{beta.est}) was seeked between \code{h.min.opt.max.mopt[1]} and \code{h.min.opt.max.mopt[3]}.}
\item{h.min.opt.max.m}{For each \eqn{\theta} in \eqn{\Theta_n}, the optimal \eqn{h} (\code{h.min.opt.max.m[m,2]}) used by the optimal \eqn{\beta} for this \eqn{\theta} was seeked between \code{h.min.opt.max.m[m,1]} and \code{h.min.opt.max.m[m,3]}.}
\item{h.seq.opt}{Sequence of eligible values for \eqn{h} considered to seek for \code{h.opt}.}
\item{theta.seq.norm}{The vector \code{theta.seq.norm[j,]} contains the coefficientes in the B-spline basis of the jth functional index in \eqn{\Theta_n}.}
\item{vn.opt}{Selected value of \code{vn}.}
\item{...}{}
}
\references{
Ait-Saidi, A., Ferraty, F., Kassa, R., and Vieu, P. (2008) Cross-validated estimations in the single-functional index model. \emph{Statistics}, \bold{42(6)}, 475--494, \doi{https://doi.org/10.1080/02331880801980377}.

Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single-index regression. \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.

Novo, S., Aneiros, G., and Vieu, P., (2021) Sparse semiparametric regression
when predictors are mixture of functional and high-dimensional variables. \emph{TEST},
\bold{30}, 481--504, \doi{https://doi.org/10.1007/s11749-020-00728-w}.

Novo, S., Aneiros, G., and Vieu, P., (2021) A kNN procedure in semiparametric
functional data analysis. \emph{Statistics and Probability Letters}, \bold{171}, 109028, \doi{https://doi.org/10.1016/j.spl.2020.109028}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
See also \code{\link{fsim.kernel.fit}}, \code{\link{predict.sfplsim.kernel}} and  \code{\link{plot.sfplsim.kernel}}

Alternative procedure \code{\link{sfplsim.kNN.fit}}.
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

#SSFPLSIM fit. Convergence errors for some theta are obtained.
ptm=proc.time()
fit<-sfplsim.kernel.fit(x=X[train,], z=z.com[train,], y=y[train],
      max.q.h=0.35,lambda.min.l=0.01,
      max.iter=5000, nknot.theta=4,criterion="BIC",nknot=20)
proc.time()-ptm

#Results
fit
names(fit)
}
}

