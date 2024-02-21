\name{PVS.fit}
\alias{PVS.fit}
\title{
PVS estimation
}
\description{
This function computes the partitioning variable selection algorithm (PVS) for sparse linear regression with covariates with functional origin.

This algorithm involves the penalised least-squares regularization procedure, which 
requires an objective criterion (\code{criterion}) to select the number of covariates in the reduced model \code{w.opt} and the penalisation parameter \code{lambda.opt}. 
}

\usage{
PVS.fit(z, y, train.1, train.2, lambda.min = NULL, lambda.min.h = NULL,
  lambda.min.l = NULL, factor.pn = 1, nlambda = 100, vn = ncol(z),
  nfolds = 10, seed = 123, wn = c(10,15,20),range.grid=NULL, criterion =
  c("GCV", "BIC", "AIC","k-fold-CV"), penalty = c("grLasso", "grMCP", 
  "grSCAD", "gel", "cMCP", "gBridge", "gLasso", "gMCP"), max.iter = 1000)
}

\arguments{
\item{z}{
Matrix containing the observations of the functional covariate that is discretised collected by row.
}
  \item{y}{
Vector containing the scalar response.
}
  \item{train.1}{
Indexes of the data used as the training sample in the 1st step. The default is \code{train.1<-1:ceiling(n/2)}.
}
  \item{train.2}{
Indexes of the data used as the training sample in the 2nd step. The default is \code{train.2<-(ceiling(n/2)+1):n}.
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
\item{range.grid}{
Vector of length 2 containing the endpoints of the grid at which the observations of the functional covariate \code{x} are evaluated (i.e. the range of the discretization). If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the size of the discretization size of \code{x} (i.e. \code{ncol(x))}.
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
The sparse linear model with covariates coming from the discretization of a curve is given by the expression
	\deqn{Y_i=\sum_{j=1}^{p_n}\beta_{0j}\zeta_i(t_j)+\varepsilon_i,\ \ \ (i=1,\dots,n)}
where 
\itemize{
\item \eqn{Y_i} is a real random response and  \eqn{\zeta_i} is supposed to be a random curve defined on some interval \eqn{[a,b]} which  is observed at the points \eqn{a\leq t_1<\dots<t_{p_n}\leq b}. 
\item  \eqn{\mathbf{\beta}_0=(\beta_{01},\dots,\beta_{0p_n})^{\top}} is a vector of unknown real coefficients.
\item \eqn{\varepsilon_i} denotes the random error.
}
In  this model, we assume that only a few scalar variables from the set \eqn{\{\zeta(t_1),\dots,\zeta(t_{p_n})\}}  form part of the model. Therefore, we must select the relevant variables (the impact points of the curve \eqn{\zeta} on the response) and estimate the model.

In this function, this model is fitted using the PVS. The PVS is an algorithm with two steps, so we split the sample into two independent subsamples (asymptotically of the same size \eqn{n_1\sim n_2\sim n/2}), one of them to be used in the first stage of the method and the other in the second stage.
\deqn{
\mathcal{E}^{\mathbf{1}}=\{(\zeta_i,\mathcal{X}_i,Y_i),\quad i=1,\dots,n_1\},
}
\deqn{
\mathcal{E}^{\mathbf{2}}=\{(\zeta_i,\mathcal{X}_i,Y_i),\quad i=n_1+1,\dots,n_1+n_2=n\}.
} 
 Note that these two subsamples are specified to the programme by means of the arguments \code{train.1} and \code{train.2}.
The superscript \eqn{\mathbf{s}} with \eqn{\mathbf{s}=\mathbf{1},\mathbf{2}} indicates the stage of the method in which the sample, function, variable or parameter is involved. 

To explain the algorithm we assume, without lost of generality, that the number \eqn{p_n} of linear covariates can be expressed as follows: \eqn{p_n=q_nw_n} with \eqn{q_n} and \eqn{w_n} integers.
\enumerate{
\item \bold{First step}. A reduced model is consider, discarding many linear covariates. The penalised least-squares procedure is applied to the reduced model using only the subsample \eqn{\mathcal{E}^{\mathbf{1}}}. Specifically:
\itemize{
\item Consider a subset of the initial \eqn{p_n} linear covariates, which contains only \eqn{w_n} equally spaced discretized observations of  \eqn{\zeta} covering the whole interval  \eqn{[a,b]}. This subset is the following:
 \deqn{
	\mathcal{R}_n^{\mathbf{1}}=\left\{\zeta\left(t_k^{\mathbf{1}}\right),\ \ k=1,\dots,w_n\right\},
} 
	where  \eqn{t_k^{\mathbf{1}}=t_{\left[(2k-1)q_n/2\right]}} and  \eqn{\left[z\right]} denotes the smallest integer not less than the real number \eqn{z}.The size (cardinal) of this subset is provided to the program in the argument \code{wn} (which contains a sequence of eligible sizes).

\item Consider the following reduced model, which involves only the \eqn{w_n} linear covariates belonging to \eqn{\mathcal{R}_n^{\mathbf{1}}}:
	\deqn{
	Y_i=\sum_{k=1}^{w_n}\beta_{0k}^{\mathbf{1}}\zeta_i(t_k^{\mathbf{1}})+\varepsilon_i^{\mathbf{1}}.
}
The penalised least-squares variable selection procedure is applied to the reduced model. This is done by means of the function \code{\link{lm.pels.fit}}, which requires the remaining arguments (for details, see the documentation of the function \code{\link{lm.pels.fit}}). The estimates obtained after that are the outputs of the first step of the algorithm.
}

\item \bold{Second step}. The variables selected in the first step and the variables in the neighbourhood of the ones selected are included. Then the penalised least-squares procedure is performed again. For that, we consider only the subsample \eqn{\mathcal{E}^{\mathbf{2}}}. Specifically:
	\itemize{
		\item Consider a new set of variables :
		\deqn{
		\mathcal{R}_n^{\mathbf{2}}=\bigcup_{\left\{k,\widehat{\beta}_{0k}^{\mathbf{1}}\not=0\right\}}\left\{\zeta(t_{(k-1)q_n+1}),\dots,\zeta(t_{kq_n})\right\}.
	}
		Denoting by \eqn{r_n=\sharp(\mathcal{R}_n^{\mathbf{2}})}, we can rename the variables in \eqn{\mathcal{R}_n^{\mathbf{2}}} as follows:
	\deqn{
		\mathcal{R}_n^{\mathbf{2}}=\left\{\zeta(t_1^{\mathbf{2}}),\dots,\zeta(t_{r_n}^{\mathbf{2}})\right\},
		}
	\item  Consider the following model, which involves only the linear covariates belonging to \eqn{\mathcal{R}_n^{\mathbf{2}}}
		\deqn{
		Y_i=\sum_{k=1}^{r_n}\beta_{0k}^{\mathbf{2}}\zeta_i(t_k^{\mathbf{2}})+\varepsilon_i^{\mathbf{2}}.}
		The penalized least-squares variable selection procedure, with kernel estimation, is applied to this model by means of the function \code{\link{lm.pels.fit}}. 
}
}
The outputs of the second step are the estimates of the model obtained with the PVS algorithm. For further details on this algorithm, see Aneiros and Vieu (2014).

\bold{Remark}: If the condition  \eqn{p_n=w_n q_n} fails, the function considers not fixed  \eqn{q_n=q_{n,k}} values \eqn{k=1,\dots,w_n},  when \eqn{p_n/w_n} is not an integer number. Specifically:
\deqn{
	q_{n,k}= \left\{\begin{array}{ll}
	[p_n/w_n]+1 &   k\in\{1,\dots,p_n-w_n[p_n/w_n]\},\\
	{[p_n/w_n]} & k\in\{p_n-w_n[p_n/w_n]+1,\dots,w_n\},
	\end{array}
	\right.
}
	where \eqn{[z]} denotes the integer part of the real number \eqn{z}.
	
}
\value{
\item{call}{The matched call.}
\item{fitted.values}{Estimated scalar response.}
\item{residuals}{Differences between \code{y} and the \code{fitted.values}}
\item{beta.est}{\eqn{\hat{\mathbf{\beta}}} (i. e. estimate of \eqn{\mathbf{\beta}_0} when the optimal tuning parameters \code{w.opt} and \code{lambda.opt} are used).}
\item{indexes.beta.nonnull}{Indexes of the non-zero \eqn{\hat{\beta_{j}}}.}
\item{w.opt}{Selected size for  \eqn{\mathcal{R}_n^{\mathbf{1}}}.}
\item{lambda.opt}{Selected value of the penalisation parameter \eqn{\lambda} (when \code{w.opt} is considered).}
\item{IC}{Value of the criterion function considered to select \code{w.opt} and \code{lambda.opt}.}
\item{beta2}{Estimate of \eqn{\mathbf{\beta}_0^{\mathbf{2}}} for each value of the sequence \code{wn}.}
\item{indexes.beta.nonnull2}{Indexes of the non-zero linear coefficients after the step 2 of the method for each value of the sequence \code{wn}.} 
\item{IC2}{Optimal value of the criterion function in the second step for each value of the sequence \code{wn}.}
\item{lambda2}{Selected value of penalisation parameter in the second step for each value of the sequence \code{wn}.}
\item{index02}{Indexes of the covariates (in the whole set of \eqn{p_n}) used to build \eqn{\mathcal{R}_n^{\mathbf{2}}} for each value of the sequence \code{wn}.} 
\item{beta1}{Estimate of \eqn{\mathbf{\beta}_0^{\mathbf{1}}} for each value of the sequence \code{wn}.}
\item{IC1}{Optimal value of the criterion function in the first step for each value of the sequence \code{wn}.}
\item{lambda1}{Selected value of penalisation parameter in the first step for each value of the sequence \code{wn}.}
\item{index01}{Indexes of the covariates (in the whole set of \eqn{p_n}) used to build \eqn{\mathcal{R}_n^{\mathbf{1}}} for each value of the sequence \code{wn}.}
\item{index1}{Indexes of the non-zero linear coefficients after the step 1 of the method for each value of the sequence \code{wn}.}
\item{...}{}
}
\references{
Aneiros, G. and Vieu, P. (2014) Variable selection in infinite-dimensional problems. \emph{Statistics  & Probability Letters}, \bold{94}, 12--20, \doi{https://doi.org/10.1016/j.spl.2014.06.025}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}



\seealso{
See also \code{\link{lm.pels.fit}}.
}
\examples{
data(Sugar)

y<-Sugar$ash
z<-Sugar$wave.240

#Outliers
index.y.25 <- y > 25
index.atip <- index.y.25
(1:268)[index.atip]


#Dataset to model
z.sug<- z[!index.atip,]
y.sug <- y[!index.atip]

train<-1:216

ptm=proc.time()
fit<- PVS.fit(z=z.sug[train,], y=y.sug[train],train.1=1:108,train.2=109:216,
        lambda.min.h=0.2,criterion="BIC", penalty="grSCAD", max.iter=5000)
proc.time()-ptm

fit 
names(fit)

 
}

