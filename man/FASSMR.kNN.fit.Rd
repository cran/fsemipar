\name{FASSMR.kNN.fit}
\alias{FASSMR.kNN.fit}
\title{
Impact point selection with FASSMR and kNN estimation
}
\description{
This function implements the Fast Algorithm for Sparse Semiparametric Multi-functional Regression (FASSMR) with kNN estimation. This algorithm is specifically designed for estimating multi-functional partial linear single-index models, which incorporate multiple scalar variables and a functional covariate as predictors. These scalar variables are derived from the discretisation of a curve and have linear effect while the functional covariate exhibits a single-index effect. 

FASSMR selects the impact points of the discretised curve and estimates the model. The algorithm employs a penalised least-squares regularisation procedure, integrated with kNN estimation using Nadaraya-Watson weights. It uses B-spline expansions to represent curves and eligible functional indexes. Additionally, it utilises an objective criterion (\code{criterion}) to determine the initial number of covariates in the reduced model (\code{w.opt}), the number of neighbours (\code{k.opt}), and the penalisation parameter (\code{lambda.opt}).
}
\usage{
FASSMR.kNN.fit(x, z, y, seed.coeff = c(-1, 0, 1), order.Bspline = 3, 
nknot.theta = 3,  knearest = NULL, min.knn = 2, max.knn = NULL, step = NULL,  
kind.of.kernel = "quad",range.grid = NULL, nknot = NULL, lambda.min = NULL, 
lambda.min.h = NULL, lambda.min.l = NULL, factor.pn = 1, nlambda = 100, 
vn = ncol(z), nfolds = 10, seed = 123, wn = c(10, 15, 20), criterion = "GCV", 
penalty = "grSCAD", max.iter = 1000, n.core = NULL)
}
\arguments{
  \item{x}{
Matrix containing the observations of the functional covariate collected by row (functional single-index component).
}
  \item{z}{
Matrix containing the observations of the functional covariate that is discretised collected by row (linear component).
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
 
    \item{kind.of.kernel}{
The type of kernel function used. Currently, only Epanechnikov kernel (\code{"quad"}) is available.
}
 \item{range.grid}{
Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretisation). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the discretisation size of \code{x} (i.e. \code{ncol(x))}.
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
  \item{vn}{
Positive integer or vector of positive integers indicating the number of groups of consecutive variables to be penalised together. The default value is \code{vn=ncol(z)}, resulting in the individual penalization of each scalar covariate.
}
  \item{nfolds}{
Positive integer indicating the number of cross-validation folds (used when \code{criterion="k-fold-CV"}). Default is 10.
}
  \item{seed}{
You may set the seed for the random number generator to ensure reproducible results (applicable when \code{criterion="k-fold-CV"} is used). The default seed value is 123.
}
  \item{wn}{
A vector of positive integers indicating the eligible number of covariates in the reduced model. For more information, refer to the section \code{Details}. The default is \code{c(10,15,20)}.
}
  \item{criterion}{
The criterion used to select the tuning and regularisation parameters: \code{wn.opt}, \code{k.opt} and \code{lambda.opt} (also \code{vn.opt} if needed). Options include \code{"GCV"}, \code{"BIC"}, \code{"AIC"}, or \code{"k-fold-CV"}. The default setting is \code{"GCV"}.
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
The multi-functional partial linear single-index model (MFPLSIM) is given by the expression
	\deqn{Y_i=\sum_{j=1}^{p_n}\beta_{0j}\zeta_i(t_j)+r\left(\left<\theta_0,X_i\right>\right)+\varepsilon_i,\ \ \ (i=1,\dots,n),}
where: 
\itemize{
\item \eqn{Y_i} is a real random response and \eqn{X_i} denotes a random element belonging to some separable Hilbert space \eqn{\mathcal{H}} with inner product denoted by \eqn{\left\langle\cdot,\cdot\right\rangle}. The second functional predictor \eqn{\zeta_i} is assumed to be a curve defined on some interval \eqn{[a,b]} which  is observed at the points \eqn{a\leq t_1<\dots<t_{p_n}\leq b}. 
\item  \eqn{\mathbf{\beta}_0=(\beta_{01},\dots,\beta_{0p_n})^{\top}} is a vector of unknown real coefficients and \eqn{r(\cdot)} denotes a smooth unknown link function. In addition, \eqn{\theta_0} is an unknown functional direction in \eqn{\mathcal{H}}.  
\item \eqn{\varepsilon_i} denotes the random error.
}
In  the MFPLSIM, we assume that only a few scalar variables from the set \eqn{\{\zeta(t_1),\dots,\zeta(t_{p_n})\}} form part of the model. Therefore, we must select the relevant variables in the linear component (the impact points of the curve \eqn{\zeta} on the response) and estimate the model.


In this function, the MFPLSIM is fitted using the FASSMR algorithm.  The main idea of this algorithm is to consider a reduced model, with only some (very few) linear covariates (but covering the entire discretization interval of \eqn{\zeta}), and discarding directly the other linear covariates (since it is expected that they contain very similar information about the response). 

To explain the algorithm, we assume, without loss of generality, that the number \eqn{p_n} of linear covariates can be expressed as follows: \eqn{p_n=q_nw_n} with \eqn{q_n} and \eqn{w_n} integers. 
This consideration allows us to build a subset of the initial \eqn{p_n} linear covariates, containging only \eqn{w_n} equally spaced discretised observations of  \eqn{\zeta} covering the entire interval  \eqn{[a,b]}. This subset is the following:
 \deqn{
	\mathcal{R}_n^{\mathbf{1}}=\left\{\zeta\left(t_k^{\mathbf{1}}\right),\ \ k=1,\dots,w_n\right\},
} 
	where  \eqn{t_k^{\mathbf{1}}=t_{\left[(2k-1)q_n/2\right]}} and  \eqn{\left[z\right]} denotes the smallest integer not less than the real number \eqn{z}.

We consider the following reduced model, which involves only the linear covariates belonging to \eqn{\mathcal{R}_n^{\mathbf{1}}}:
	\deqn{
	Y_i=\sum_{k=1}^{w_n}\beta_{0k}^{\mathbf{1}}\zeta_i(t_k^{\mathbf{1}})+r^{\mathbf{1}}\left(\left<\theta_0^{\mathbf{1}},\mathcal{X}_i\right>\right)+\varepsilon_i^{\mathbf{1}}.
}
The program receives the eligible numbers of linear covariates for building the reduced model through the argument \code{wn}.
Then, the penalised least-squares variable selection procedure, with kNN estimation, is applied to the reduced model. This is done using the function \code{\link{sfplsim.kNN.fit}}, which requires the remaining arguments (for details, see the documentation of the function \code{\link{sfplsim.kNN.fit}}). The estimates obtained are the outputs of the FASSMR algorithm. For further details on this algorithm, see Novo et al. (2021).

\bold{Remark}: If the condition  \eqn{p_n=w_n q_n} is not met (then \eqn{p_n/w_n} is not an integer number), the function considers variable  \eqn{q_n=q_{n,k}} values \eqn{k=1,\dots,w_n}. Specifically:
\deqn{
	q_{n,k}= \left\{\begin{array}{ll}
	[p_n/w_n]+1 &   k\in\{1,\dots,p_n-w_n[p_n/w_n]\},\\
	{[p_n/w_n]} & k\in\{p_n-w_n[p_n/w_n]+1,\dots,w_n\},
	\end{array}
	\right.
}
	where \eqn{[z]} denotes the integer part of the real number \eqn{z}.

The function supports parallel computation. To avoid it, we can set \code{n.core=1}.
}
\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}.}
\item{beta.est}{\eqn{\hat{\mathbf{\beta}}} (i.e. estimate of \eqn{\mathbf{\beta}_0} when the optimal tuning parameters \code{w.opt}, \code{lambda.opt}, \code{k.opt} and \code{vn.opt} are used).}
\item{beta.red}{Estimate of \eqn{\beta_0^{\mathbf{1}}} in the reduced model when the optimal tuning parameters \code{w.opt}, \code{lambda.opt}, \code{k.opt} and \code{vn.opt} are used.}
\item{theta.est}{Coefficients of \eqn{\hat{\theta}} in the B-spline basis (i.e. estimate of \eqn{\theta_0} when the optimal tuning parameters \code{w.opt}, \code{lambda.opt}, \code{k.opt} and \code{vn.opt} are used): a vector of \code{length(order.Bspline+nknot.theta)}.}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{k.opt}{Selected number of nearest neighbours (when \code{w.opt} is considered).}
\item{w.opt}{Selected size for  \eqn{\mathcal{R}_n^{\mathbf{1}}}.}
\item{lambda.opt}{Selected value for the penalisation parameter (when \code{w.opt} is considered).}
\item{IC}{Value of the criterion function considered to select \code{w.opt}, \code{lambda.opt}, \code{k.opt} and \code{vn.opt}.}
\item{vn.opt}{Selected value of \code{vn} (when \code{w.opt} is considered).}
\item{beta.w}{Estimate of \eqn{\beta_0^{\mathbf{1}}} for each value of the sequence \code{wn} (i.e. for each number of covariates in the reduced model).}
\item{theta.w}{Estimate of \eqn{\theta_0^{\mathbf{1}}} for each value of the sequence \code{wn} (i.e. its coefficients in the B-spline basis).} 
\item{IC.w}{Value of the criterion function for each value of the sequence \code{wn}.}
\item{indexes.beta.nonnull.w}{Indexes of the non-zero linear coefficients for each value of the sequence \code{wn}.} 
\item{lambda.w}{Selected value of penalisation parameter for each value of the sequence \code{wn}.}
\item{k.w}{Selected number of neighbours for each value of the sequence \code{wn}.}
\item{index01}{Indexes of the covariates (in the entire set of \eqn{p_n}) used to build  \eqn{\mathcal{R}_n^{\mathbf{1}}} for each value of the sequence \code{wn}.}
\item{...}{}
}
\references{
Novo, S., Vieu, P., and Aneiros, G., (2021) Fast and efficient algorithms for
sparse semiparametric bi-functional regression. \emph{Australian and New Zealand
Journal of Statistics}, \bold{63}, 606--638, \doi{https://doi.org/10.1111/anzs.12355}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
See also \code{\link{sfplsim.kNN.fit}, \link{predict.FASSMR.kNN}}, \code{\link{plot.FASSMR.kNN}} and \code{\link{IASSMR.kNN.fit}}.

Alternative method \code{\link{FASSMR.kernel.fit}}
}
\examples{
\donttest{
data(Sugar)


y<-Sugar$ash
x<-Sugar$wave.290
z<-Sugar$wave.240

#Outliers
index.y.25 <- y > 25
index.atip <- index.y.25
(1:268)[index.atip]


#Dataset to model
x.sug <- x[!index.atip,]
z.sug<- z[!index.atip,]
y.sug <- y[!index.atip]

train<-1:216
ptm=proc.time()
fit<- FASSMR.kNN.fit(x=x.sug[train,],z=z.sug[train,], y=y.sug[train], 
        nknot.theta=2, lambda.min.l=0.03, max.knn=20,nknot=20,criterion="BIC",
        max.iter=5000)
proc.time()-ptm

fit
names(fit)
}
  
}  
