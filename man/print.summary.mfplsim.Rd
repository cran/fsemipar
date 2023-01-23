\name{print.summary.mfplsim}
\alias{print.FASSMR.kernel}
\alias{print.FASSMR.kNN}
\alias{print.IASSMR.kernel}
\alias{print.IASSMR.kNN}
\alias{summary.FASSMR.kernel}
\alias{summary.FASSMR.kNN}
\alias{summary.IASSMR.kernel}
\alias{summary.IASSMR.kNN}
\title{
Summarize information from multi-functional partial linear single-index model (MFPLSIM) estimation methods
}
\description{
\code{summary} and \code{print} functions for \code{FASSMR.kernel.fit}, \code{FASSMR.kNN.fit}, \code{IASSMR.kernel.fit} and \code{IASSMR.kNN.fit}.
}
\usage{
\method{print}{FASSMR.kernel}(x, \dots)
\method{print}{FASSMR.kNN}(x, \dots)
\method{print}{IASSMR.kernel}(x, \dots)
\method{print}{IASSMR.kNN}(x, \dots)
\method{summary}{FASSMR.kernel}(object, \dots)
\method{summary}{FASSMR.kNN}(object, \dots)
\method{summary}{IASSMR.kernel}(object, \dots)
\method{summary}{IASSMR.kNN}(object, \dots)
}
\arguments{
  \item{x}{
Output of the \code{FASSMR.kernel.fit}, \code{FASSMR.kNN.fit}, \code{IASSMR.kernel.fit} or \code{IASSMR.kNN.fit} functions (i.e. an object of the class \code{FASSMR.kernel},  \code{FASSMR.kNN}, \code{IASSMR.kernel} or \code{IASSMR.kNN}).
}
  \item{\dots}{
}
\item{object}{
Output of the \code{FASSMR.kernel.fit}, \code{FASSMR.kNN.fit}, \code{IASSMR.kernel.fit} or \code{IASSMR.kNN.fit} functions (i.e. an object of the class \code{FASSMR.kernel},  \code{FASSMR.kNN}, \code{IASSMR.kernel} or \code{IASSMR.kNN}).
}
}

\value{
\itemize{
\item{The matched call.}
\item{The optimal value of the tunning parameter (\code{h.opt} or \code{k.opt}).}
\item{The optimal initial number of covariates to build the reduced model (\code{w.opt}).}
\item{Coefficients of \eqn{\hat{\theta}} in the B-spline basis (\code{theta.est}): a vector of \code{length(order.Bspline+nknot.theta).}}
\item{The estimated vector of linear coefficients (\code{beta.est}).}
\item{The number of non-zero components in \code{beta.est}.}
\item{The indexes of the non-zero components in \code{beta.est}.}
\item{The optimal value of the penalisation parameter (\code{lambda.opt}).}
\item{The optimal value of the criterion function, i. e. the value  obtained with \code{w.opt}, \code{lambda.opt}, \code{vn.opt}  and \code{h.opt}/\code{k.opt}}
\item{Minimum value of the penalized profile least-squares function. That is, the value obtained using \code{theta.est} and \code{beta.est}.}
\item{The penalty function used.}
\item{The criterion used to select the initial number of covariates used to build the reduced model, the tunning parameter, the penalisation parameter and \code{Vn}.}


}
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
\code{FASSMR.kernel.fit}, \code{FASSMR.kNN.fit}, \code{IASSMR.kernel.fit} and \code{IASSMR.kNN.fit}.
}
