\name{lm.pels.fit}
\alias{lm.pels.fit}

\title{
Regularised fit of sparse linear regression
}
\description{
This function fits a sparse linear model between a scalar response and a vector of scalar covariates. It employs a penalised least-squares regularisation procedure, with either (group)SCAD or (group)LASSO penalties. The method utilises an objective criterion (\code{criterion}) to select the optimal regularisation parameter (\code{lambda.opt}).
}
\usage{
lm.pels.fit(z, y, lambda.min = NULL, lambda.min.h = NULL, lambda.min.l = NULL,
factor.pn = 1, nlambda = 100, lambda.seq = NULL, vn = ncol(z), nfolds = 10, 
seed = 123, criterion = "GCV", penalty = "grSCAD", max.iter = 1000)
}

\arguments{
   \item{z}{
Matrix containing the observations of the covariates collected by row.
}
  \item{y}{
Vector containing the scalar response.
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
The criterion used to select the regularisation parameter \code{lambda.opt} (also \code{vn.opt} if needed). Options include \code{"GCV"}, \code{"BIC"}, \code{"AIC"}, or \code{"k-fold-CV"}. The default setting is \code{"GCV"}.
}
  \item{penalty}{
The penalty function applied in the penalised least-squares procedure. Currently, only "grLasso" and "grSCAD" are implemented. The default is "grSCAD".
}
 
  \item{max.iter}{
Maximum number of iterations allowed across the entire path. The default value is 1000.
}
}
\details{
The sparse linear model (SLM) is given by the expression:
\deqn{
Y_i=Z_{i1}\beta_{01}+\dots+Z_{ip_n}\beta_{0p_n}+\varepsilon_i\ \ \ i=1,\dots,n,
}
where \eqn{Y_i} denotes a scalar response, \eqn{Z_{i1},\dots,Z_{ip_n}} are real covariates. In this equation, \eqn{\mathbf{\beta}_0=(\beta_{01},\dots,\beta_{0p_n})^{\top}} is a vector of unknown real parameters and \eqn{\varepsilon_i} represents the random error.

In this function, the SLM is fitted using a penalised least-squares (PeLS) approach by minimising 
\deqn{
\mathcal{Q}\left(\mathbf{\beta}\right)=\frac{1}{2}\left(\mathbf{Y}-\mathbf{Z}\mathbf{\beta}\right)^{\top}\left(\mathbf{Y}-\mathbf{Z}\mathbf{\beta}\right)+n\sum_{j=1}^{p_n}\mathcal{P}_{\lambda_{j_n}}\left(|\beta_j|\right), \quad (1)
}
where \eqn{\mathbf{\beta}=(\beta_1,\ldots,\beta_{p_n})^{\top}, \ \mathcal{P}_{\lambda_{j_n}}\left(\cdot\right)} is a penalty function (specified in the argument \code{penalty}) and \eqn{\lambda_{j_n} > 0} is a tuning parameter.
To reduce the number of tuning parameters, \eqn{\lambda_j}, to be selected for each sample, we consider \eqn{\lambda_j = \lambda \widehat{\sigma}_{\beta_{0,j,OLS}}}, where \eqn{\beta_{0,j,OLS}} denotes the OLS estimate of \eqn{\beta_{0,j}} and \eqn{\widehat{\sigma}_{\beta_{0,j,OLS}}} is the estimated standard deviation. The parameter \eqn{\lambda} is selected using the objetive criterion specified in the argument \code{criterion}.

For further details on the estimation procedure of the SLM, see e.g. Fan and Li. (2001). The PeLS objective function  is minimised using the R function \code{grpreg} of the package \code{grpreg} (Breheny and Huang, 2015).

\bold{Remark}:  It should be noted that if  we set \code{lambda.seq} to \eqn{=0}, we obtain the non-penalised estimation of the model, i.e. the OLS estimation. Using \code{lambda.seq} with a vaule \eqn{\not=0} is advisable when suspecting the presence of irrelevant variables.
}

\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}.}
\item{beta.est}{Estimate of \eqn{\beta_0} when the optimal penalisation parameter \code{lambda.opt} and \code{vn.opt} are used.}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{lambda.opt}{Selected value of lambda.}
\item{IC}{Value of the criterion function considered to select \code{lambda.opt} and \code{vn.opt}.}
\item{vn.opt}{Selected value of \code{vn}.}
\item{...}{}
}
\references{
Breheny, P., and Huang, J. (2015) Group descent algorithms for nonconvex penalized linear and
logistic regression models with grouped predictors. \emph{Statistics and Computing}, \bold{25}, 173--187, \doi{https://doi.org/10.1007/s11222-013-9424-2}.

Fan, J., and Li, R. (2001) Variable selection via nonconcave penalized
likelihood and its oracle properties. \emph{Journal of the American Statistical Association}, \bold{96}, 1348--1360, \doi{https://doi.org/10.1198/016214501753382273}.
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

#LM fit 
ptm=proc.time()
fit<-lm.pels.fit(z=z.com[train,], y=y[train],lambda.min.h=0.02,
      lambda.min.l=0.01,factor.pn=2, max.iter=5000, criterion="BIC")
proc.time()-ptm

#Results
fit
names(fit)
 
}

