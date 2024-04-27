\name{print.summary.fsim}
\alias{print.fsim.kernel}
\alias{print.fsim.kNN}
\alias{summary.fsim.kernel}
\alias{summary.fsim.kNN}
\title{
Summarise information from FSIM estimation 
}
\description{
\code{summary} and \code{print} functions for \code{fsim.kNN.fit}, \code{fsim.kNN.fit.optim}, \code{fsim.kernel.fit} and \code{fsim.kernel.fit.optim}.
}
\usage{

\method{print}{fsim.kernel}(x, \dots)
\method{print}{fsim.kNN}(x, \dots)
\method{summary}{fsim.kernel}(object, \dots)
\method{summary}{fsim.kNN}(object, \dots)
}
\arguments{
  \item{x}{
Output of the \code{fsim.kernel.fit}, \code{fsim.kernel.fit.optim}, \code{fsim.kNN.fit} or \code{fsim.kNN.fit.optim} functions (i.e. an object of the class \code{fsim.kernel} or \code{fsim.kNN}).
}
\item{\dots}{ Further arguments.}
\item{object}{
Output of the \code{fsim.kernel.fit}, \code{fsim.kernel.fit.optim}, \code{fsim.kNN.fit} or \code{fsim.kNN.fit.optim} functions (i.e. an object of the class \code{fsim.kernel} or \code{fsim.kNN}).
}
}

\value{
\itemize{
\item{The matched call.}
\item{The optimal value of the tunning parameter (\code{h.opt} or \code{k.opt}).}
\item{Coefficients of \eqn{\hat{\theta}} in the B-spline basis (\code{theta.est}: a vector of \code{length(order.Bspline+nknot.theta).}}
\item{Minimum value of the CV function, i.e. the value of CV for \code{theta.est} and \code{h.opt}/\code{k.opt}.}
\item{R squared.}
\item{Residual variance.}
\item{Residual degrees of freedom.}

}
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}

\seealso{
\code{fsim.kernel.fit} and \code{fsim.kNN.fit}.
}
