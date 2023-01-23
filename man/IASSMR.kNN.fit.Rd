\name{IASSMR.kNN.fit}
\alias{IASSMR.kNN.fit}
\title{
IASSMR with kNN estimation
}
\description{
This function computes the improved algorithm for sparse semiparametric multi-functional regression (IASSMR) with kNN estimation.

This algorithm involves the penalised least-squares regularization procedure combined with \eqn{k}-nearest neighbours (kNN) estimation with Nadaraya-Watson weights.
The procedure requires the  B-spline representation to estimate the functional index \eqn{\theta_0} and an objective criterion (\code{criterion}) to select the number of covariates in the reduced model (\code{w.opt}), the number of neighbours (\code{k.opt}) and the penalisation parameter (\code{lambda.opt}). 
}
\usage{
IASSMR.kNN.fit(x, z, y, train.1=NULL, train.2=NULL, seed.coeff = c(-1, 0, 1), 
  order.Bspline = 3, nknot.theta = 3, t0 = NULL, knearest = NULL, 
  min.knn = 2, max.knn = NULL, step = NULL, range.grid = NULL, 
  kind.of.kernel = "quad", nknot = NULL, lambda.min = NULL, 
  lambda.min.h = NULL, lambda.min.l = NULL, factor.pn = 1,
  nlambda = 100,vn = ncol(z), nfolds = 10, seed = 123, wn = c(10, 15, 20),
  criterion = c("GCV", "BIC", "AIC", "k-fold-CV"), 
  penalty = c("grLasso", "grMCP", "grSCAD", "gel", "cMCP",
  "gBridge", "gLasso", "gMCP"), max.iter = 1000)
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
  \item{train.1}{Indexes of the data used as the training sample in the 1st step. The default is  \code{train.1<-1:ceiling(n/2)}.
}
  \item{train.2}{Indexes of the data used as the training sample in the 2nd step.  The default is \code{train.2<-(ceiling(n/2)+1):n}.
}
  \item{seed.coeff}{
Vector of initial values used to  build the set \eqn{\Theta_n} (see section \code{Details}). The coefficients for the B-spline representation of each eligible functional index \eqn{\theta \in \Theta_n} are obtained from \code{seed.coeff}.  The default is \code{c(-1,0,1)}.
}
  \item{order.Bspline}{
Positive integer giving the order of the B-spline basis functions. This is the number of coefficients in each piecewise polynomial segment. The default is 3.
}
  \item{nknot.theta}{
Positive integer indicating the number of uniform interior knots of the B-spline basis for the B-spline representation of \eqn{\theta_0}. The default is 3.
}
  \item{t0}{
Value in the domain of the functional indexes at which we evaluate them to  build the set \eqn{\Theta_n}. We assume \eqn{\theta_0(t_0)>0} for some arbitrary \eqn{t_0} in the domain to ensure model identifiability. If \code{t0=NULL}, then \code{mean(range.grid)} is considered.
}
  \item{knearest}{
Vector of positive integers containing the sequence in which the  number of nearest neighbours \code{k.opt} is selected. If \code{knearest=NULL}, then \code{knearest <- seq(from =min.knn, to = max.knn, by = step)}.
}
  \item{min.knn}{
Positive integer indicating the minumum value of the sequence in which the  number of nearest neighbours \code{k.opt} is selected (thus, this number must be smaller than the sample size). The default is 2.
}
  \item{max.knn}{
Positive integer indicating the maximum value of the sequence in which the number of nearest neighbours \code{k.opt} is selected (thus, this number must be larger than \code{min.kNN} and smaller than the sample size). The default is \code{max.knn <- n\%/\%2}, being \eqn{n=n_1} in the 1st step and \eqn{n=n_2} in the 2nd step of the method (see section \code{Details}).
}
  \item{step}{
Positive integer used to build the sequence of k-nearest neighbours in the following way: \code{min.knn, min.knn + step, min.knn + 2*step, min.knn + 3*step,...} The default is \code{step<-ceiling(n/100)}.
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
 Positive integer used to set \code{lambda.min}. The default value is 1.
}
  \item{nlambda}{
Positive integer indicating the number of values of the sequence in which \code{lambda.opt} is selected. The default is 100.
}
  \item{vn}{Positive integer or vector of positive integers indicating the number of groups of consecutive variables to be penalised together. The default value is \code{vn=ncol(z)}, which leads to the individual penalisation of each scalar covariate.
}
  \item{nfolds}{
Positive integer indicating the number of cross-validation folds (used if \code{criterion="k-fold-CV"}). Default is 10.
}
  \item{seed}{
You may set the seed of the random number generator to obtain reproducible results (used if \code{criterion="k-fold-CV"}). Default is 123.
}
  \item{wn}{
A vector of positive integers indicating the eligible number of covariates of the reduced model. See the section \code{Details}. The default is \code{c(10,15,20)}.
}
  \item{criterion}{
The criterion by which to select the regularization parameter \code{lambda.opt} and \code{k.opt}. One of \code{"GCV", "BIC", "AIC"} or \code{"k-fold-CV"}. The default is \code{"GCV"}.
}
  \item{penalty}{
The penalty function to be applied in the penalized least squares procedure. Only "grLasso" and "grSCAD" are implemented.
}
 
  \item{max.iter}{
Maximum number of iterations (total across entire path). Default is 1000.
}
}

\details{
The multi-functional partial linear single-index model (MFPLSIM) is given by the expression
	\deqn{Y_i=\sum_{j=1}^{p_n}\beta_{0j}\zeta_i(t_j)+r\left(\left<\theta_0,X_i\right>\right)+\varepsilon_i,\ \ \ (i=1,\dots,n)}
where 
\itemize{
\item \eqn{Y_i} is a real random response and \eqn{X_i} denotes a random element belonging to some separable Hilbert space \eqn{\mathcal{H}} with inner product denoted by \eqn{\left\langle\cdot,\cdot\right\rangle}. The second functional predictor \eqn{\zeta_i} is supposed to be a random curve defined on some interval \eqn{[a,b]} which  is observed at the points \eqn{a\leq t_1<\dots<t_{p_n}\leq b}. 
\item  \eqn{\pmb{\beta}_0=(\beta_{01},\dots,\beta_{0p_n})^{\top}} is a vector of unknown real coefficients and \eqn{r(\cdot)} denotes a smooth unknown link function. In addition, \eqn{\theta_0} is an unknown functional index in \eqn{\mathcal{H}}.  
\item \eqn{\varepsilon_i} denotes the random error.
}
In  the MFPLSIM, we assume that only a few scalar variables from the set \eqn{\{\zeta(t_1),\dots,\zeta(t_{p_n})\}}  form part of the model. Therefore, we must select the relevant variables in the linear component  (the impact points of the curve \eqn{\zeta} on the response) and estimate the model.

In this function, the MFPLSIM is fitted using the IASSMR. The IASSMR  (version of the PVS algorithm for semiparametric regression) is an algorithm with two steps, so we split the sample into two independent subsamples (asymptotically of the same size \eqn{n_1\sim n_2\sim n/2}), one of them to be used in the first stage of the method and the other in the second stage.
\deqn{
\mathcal{E}^{\pmb{1}}=\{(\zeta_i,\mathcal{X}_i,Y_i),\quad i=1,\dots,n_1\},
}
\deqn{
\mathcal{E}^{\pmb{2}}=\{(\zeta_i,\mathcal{X}_i,Y_i),\quad i=n_1+1,\dots,n_1+n_2=n\}.
} 
 Note that these two subsamples are specified to the programme by means of the arguments \code{train.1} and \code{train.2}.
The superscript \eqn{\pmb{s}} with \eqn{\pmb{s}=\pmb{1},\pmb{2}} indicates the stage of the method in which the sample, function, variable or parameter is involved. 

To explain the algorithm we assume, without lost of generality, that the number \eqn{p_n} of linear covariates can be expressed as follows: \eqn{p_n=q_nw_n} with \eqn{q_n} and \eqn{w_n} integers.
\enumerate{
\item \bold{First step}. The fast algorithm for sparse semiparametric multi-functional regression (FASSMR) combined with kNN estimation is applied using only the subsample \eqn{\mathcal{E}^{\pmb{1}}} (see the documentation of the function \code{\link{FASSMR.kNN.fit}}). Specifically:
\itemize{
\item Consider a subset of the initial \eqn{p_n} linear covariates, which contains only \eqn{w_n} equally spaced discretized observations of  \eqn{\zeta} covering the whole interval  \eqn{[a,b]}. This subset is the following:
 \deqn{
	\mathcal{R}_n^{\pmb{1}}=\left\{\zeta\left(t_k^{\pmb{1}}\right),\ \ k=1,\dots,w_n\right\},
} 
	where  \eqn{t_k^{\pmb{1}}=t_{\left[(2k-1)q_n/2\right]}} and  \eqn{\left[z\right]} denotes the smallest integer not less than  \eqn{z\in\mathbb{R}}.The size (cardinal) of this subset is provided to the program in the argument \code{wn} (which contains a sequence of eligible sizes).

\item Consider the following reduced model, which involves only the \eqn{w_n} linear covariates belonging to \eqn{\mathcal{R}_n^{\pmb{1}}}:
	\deqn{
	Y_i=\sum_{k=1}^{w_n}\beta_{0k}^{\pmb{1}}\zeta_i(t_k^{\pmb{1}})+r^{\pmb{1}}\left(\left<\theta_0^{\pmb{1}},X_i\right>\right)+\varepsilon_i^{\pmb{1}}.
}
The penalised least-squares variable selection procedure, with kNN estimation, is applied to the reduced model. This is done by means of the function \code{\link{sfplsim.kNN.fit}}, which requires the remaining arguments (for details, see the documentation of the function \code{\link{sfplsim.kNN.fit}}). The estimates obtained after that are the outputs of the first step of the algorithm.
}

\item \bold{Second step}. The variables selected in the first step and the variables in the neighbourhood of the ones selected are included. Then the penalised least-squares procedure, combined with kNN estimation, is carried out again. For that, we consider only the subsample \eqn{\mathcal{E}^{\pmb{2}}}. Specifically:
	\itemize{
		\item Consider a new set of variables:
		\deqn{
		\mathcal{R}_n^{\pmb{2}}=\bigcup_{\left\{k,\widehat{\beta}_{0k}^{\pmb{1}}\not=0\right\}}\left\{\zeta(t_{(k-1)q_n+1}),\dots,\zeta(t_{kq_n})\right\}.
	}
		Denoting by \eqn{r_n=\sharp(\mathcal{R}_n^{\pmb{2}})}, we can rename the variables in \eqn{\mathcal{R}_n^{\pmb{2}}} as follows:
	\deqn{
		\mathcal{R}_n^{\pmb{2}}=\left\{\zeta(t_1^{\pmb{2}}),\dots,\zeta(t_{r_n}^{\pmb{2}})\right\},
		}
	\item  Consider the following model, which involves only the linear covariates belonging to \eqn{\mathcal{R}_n^{\pmb{2}}}
		\deqn{
		Y_i=\sum_{k=1}^{r_n}\beta_{0k}^{\pmb{2}}\zeta_i(t_k^{\pmb{2}})+r^{\pmb{2}}\left(\left<\theta_0^{\pmb{2}},X_i\right>\right)+\varepsilon_i^{\pmb{2}}.}
		The penalized least-squares variable selection procedure, with kNN estimation, is applied to this model by means of the function \code{\link{sfplsim.kNN.fit}}. 
}
}
The outputs of the second step are the estimates of the MFPLSIM obtained with the IASSMR algorithm. For further details on this algorithm, see Novo et al. (2021).

\bold{Remark}: If the condition  \eqn{p_n=w_n q_n} fails, the function considers not fixed  \eqn{q_n=q_{n,k}} values \eqn{k=1,\dots,w_n},  when \eqn{p_n/w_n} is not an integer number. Specifically:
\deqn{
	q_{n,k}= \left\{\begin{aligned}
	&[p_n/w_n]+1&  \quad k\in\{1,\dots,p_n-w_n[p_n/w_n]\},\\
	&[p_n/w_n]&\quad k\in\{p_n-w_n[p_n/w_n]+1,\dots,w_n\}, \\
	\end{aligned}
	\right.
}
	where \eqn{[z]} denotes the integer part of \eqn{z\in\mathbb{R}}.
	


}


\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}}
\item{beta.est}{\eqn{\hat{\pmb{\beta}}} (i.e. estimate of \eqn{\pmb{\beta}_0} when the optimal tuning parameters \code{w.opt}, \code{lambda.opt}, \code{vn.opt} and \code{k.opt} are used).}
\item{theta.est}{Coefficients of \eqn{\hat{\theta}} in the B-spline basis (i. e. estimate of \eqn{\theta_0}when the optimal tuning parameters \code{w.opt}, \code{lambda.opt}, \code{vn.opt} and \code{k.opt} are used): a vector of \code{length(order.Bspline+nknot.theta)}.}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{k.opt}{Selected number of nearest neighbours (when \code{w.opt} is considered).}
\item{w.opt}{Selected initial number of covariates in the reduced model.}
\item{lambda.opt}{Selected value of the penalisation parameter \eqn{\lambda} (when \code{w.opt} is considered).}
\item{IC}{Value of the criterion function considered to select \code{w.opt}, \code{lambda.opt}, \code{vn.opt} and \code{k.opt}.}
\item{vn.opt}{Selected value of \code{vn} in the second step (when \code{w.opt} is considered).}
\item{beta2}{Estimate of \eqn{\pmb{\beta}_0^{\pmb{2}}} for each value of the sequence \code{wn}.}
\item{theta2}{Estimate of \eqn{\theta_0^{\pmb{2}}} for each value of the sequence \code{wn} (i.e. its coefficients in the B-spline basis).} 
\item{indexes.beta.nonnull2}{Indexes of the non-zero linear coefficients after the step 2 of the method for each value of the sequence \code{wn}.} 
\item{knn2}{Selected number of neighbours in the second step of the algorithm for each value of the sequence \code{wn}.} 
\item{IC2}{Optimal value of the criterion function in the second step for each value of the sequence \code{wn}.}
\item{lambda2}{Selected value of penalisation parameter in the second step for each value of the sequence \code{wn}.}
\item{index02}{Indexes of the covariates (in the whole set of \eqn{p_n}) used to build \eqn{\mathcal{R}_n^{\pmb{2}}} for each value of the sequence \code{wn}.} 
\item{beta1}{Estimate of \eqn{\pmb{\beta}_0^{\pmb{1}}} for each value of the sequence \code{wn}.}
\item{theta1}{Estimate of \eqn{\theta_0^{\pmb{1}}} for each value of the sequence \code{wn} (i.e. its coefficients in the B-spline basis).}
\item{knn1}{Selected number of neighbours in the first step of the algorithm for each value of the sequence \code{wn}.}
\item{IC1}{Optimal value of the criterion function in the first step for each value of the sequence \code{wn}.}
\item{lambda1}{Selected value of penalisation parameter in the first step for each value of the sequence \code{wn}.}
\item{index01}{Indexes of the covariates (in the whole set of \eqn{p_n}) used to build \eqn{\mathcal{R}_n^{\pmb{1}}} for each value of the sequence \code{wn}.}
\item{index1}{Indexes of the non-zero linear coefficients after the step 1 of the method for each value of the sequence \code{wn}.}
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
See also \code{\link{sfplsim.kNN.fit}, \link{predict.IASSMR.kNN}}, \code{\link{plot.IASSMR.kNN}} and \code{\link{FASSMR.kNN.fit}}.

Alternative method \code{\link{IASSMR.kernel.fit}}
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
fit<- IASSMR.kNN.fit(x=x.sug[train,],z=z.sug[train,], y=y.sug[train],
        train.1=1:108,train.2=109:216,nknot.theta=2,lambda.min.h=0.07, 
        lambda.min.l=0.07, max.knn=20, nknot=20,criterion="BIC", 
        penalty="grSCAD", max.iter=5000)
proc.time()-ptm

fit 
names(fit)
}
    
 
}
