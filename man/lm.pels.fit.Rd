\name{lm.pels.fit}
\alias{lm.pels.fit}
\title{
Linear model fit
}
\description{
This function fits a sparse  linear model between a scalar response and a vector of scalar covariates. 
The function uses the penalised least-squares regularization procedure.
The method requires an objective criterion (\code{criterion}) to select the regularization parameter (\code{lambda.opt}).
}
\usage{
lm.pels.fit(z, y, lambda.min = NULL, lambda.min.h = NULL, 
  lambda.min.l = NULL, factor.pn =1, nlambda = 100,lambda.seq = NULL,
  vn = ncol(z), nfolds = 10, seed = 123, criterion = c("GCV", "BIC", 
  "AIC", "k-fold-CV"), penalty = c("grLasso", "grMCP",
  "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"), 
  max.iter = 1000)
}

\arguments{
  \item{z}{
Matrix containing the observations of the covariates collected by row.
}
  \item{y}{
Vector containing the scalar response.
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
The criterion by which to select the regularization parameter \code{lambda.opt} and \code{k.opt}. One of \code{"GCV", "BIC", "AIC"} or \code{"k-fold-CV"}. The default is \code{"GCV"}.
}
  \item{penalty}{
The penalty function to be applied in the penalized least squares procedure. Only "grLasso" and "grSCAD" are implemented. 
}
  
 \item{max.iter}{
 Maximum number of iterations (total across entire path). The default is 1000.
}

}
\details{
The sparse linear model (SLM) is given by the expression:
\deqn{
Y_i=Z_{i1}\beta_{01}+\dots+Z_{ip_n}\beta_{0p_n}+\varepsilon_i\ \ \ i=1,\dots,n,
}
where \eqn{Y_i} denotes a scalar response, \eqn{Z_{i1},\dots,Z_{ip_n}} are random covariates taking values in \eqn{\mathbb{R}}. In this equation, \eqn{\pmb{\beta}_0=(\beta_{01},\dots,\beta_{0p_n})^{\top} \in \mathbb{R}^{p_n}} is a vector of unknown real parameters and \eqn{\varepsilon_i} is the random error.

In this function, the SLM is fitted using the penalised least-squares approach by minimising 
\deqn{
\mathcal{Q}\left(\pmb{\beta}\right)=\frac{1}{2}\left(\pmb{Y}-\pmb{Z}\pmb{\beta}\right)^{\top}\left(\pmb{Y}-\pmb{Z}\pmb{\beta}\right)+n\sum_{j=1}^{p_n}\mathcal{P}_{\lambda_{j_n}}\left(|\beta_j|\right), \quad (1)
}

where \eqn{\pmb{\beta}=(\beta_1,\ldots,\beta_{p_n})^{\top}, \ \mathcal{P}_{\lambda_{j_n}}\left(\cdot\right)} is a penalty function (specified in the argument \code{penalty}) and \eqn{\lambda_{j_n} > 0} is a tuning parameter.
To reduce  the quantity of tuning parameters, \eqn{\lambda_j}, to be selected for each sample, we consider \eqn{\lambda_j = \lambda \widehat{\sigma}_{\beta_{0,j,OLS}}}, where \eqn{\beta_{0,j,OLS}} denotes the OLS estimate of \eqn{\beta_{0,j}} and \eqn{\widehat{\sigma}_{\beta_{0,j,OLS}}} is the estimated standard deviation; \eqn{\lambda} is selected using the objetive criterion specified in the argument \code{criterion}.

For further details on the estimation procedure of the SLM, see, for instance, Fan and Li. (2001) or Fan and Lv (2011).

\bold{Remark}: We should note that if we set \code{lambda.seq}\eqn{=0}, we can obtain the non-penalised estimation of the model, i.e. the OLS estimation. It is convenient to use \code{lambda.seq}\eqn{\not=0} when one suspects there are irrelevant variables.
}
\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}}
\item{beta.est}{Estimate of \eqn{\beta_0} when the optimal penalisation parameter \code{lambda.opt} and \code{vn.opt} are used.}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{lambda.opt}{Selected value of lambda.}
\item{IC}{Value of the criterion function considered to select \code{lambda.opt} and \code{vn.opt}.}
\item{vn.opt}{Selected value of \code{vn}.}
\item{...}{}
}
\references{
Fan, J., and Li, R. (2001) Variable selection via nonconcave penalized
likelihood and its oracle properties. \emph{Journal of the American Statistical Association}, \bold{96}, 1348--1360, \doi{ https://doi.org/10.1198/016214501753382273}.


Fan, J., and  Lv, J. (2011) Nonconcave penalized likelihood with NP-dimensionality, \emph{IEEE Transactions on Information Theory}, \bold{57(8)}, 5467--5484, \url{https://ieeexplore.ieee.org/document/5961830}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
See also \code{\link{PVS.fit}}.
}
\examples{

data("Tecator")
y<-Tecator$fat
z1<-Tecator$protein       
z2<-Tecator$moisture

#Quadratic, cubic and interaction effects of the scalar covariates.
z.com<-cbind(z1,z2,z1^2,z2^2,z1^3,z2^3,z1*z2)
train<-1:160


#LM fit. 
ptm=proc.time()
fit<-lm.pels.fit(z=z.com[train,], y=y[train],lambda.min.h=0.02,
      lambda.min.l=0.01,factor.pn=2, max.iter=5000, criterion="BIC",
      penalty="grSCAD")
proc.time()-ptm

#Results
fit
names(fit)


  
}
