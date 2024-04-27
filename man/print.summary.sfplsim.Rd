\name{print.summary.sfplsim}
\alias{print.sfplsim.kernel}
\alias{print.sfplsim.kNN}
\alias{summary.sfplsim.kernel}
\alias{summary.sfplsim.kNN}
\title{
Summarise information from SFPLSIM estimation 
}
\description{
\code{summary} and \code{print} functions for \code{sfplsim.kNN.fit} and \code{sfplsim.kernel.fit}.
}
\usage{
\method{print}{sfplsim.kernel}(x, \dots)
\method{print}{sfplsim.kNN}(x, \dots)
\method{summary}{sfplsim.kernel}(object, \dots)
\method{summary}{sfplsim.kNN}(object, \dots)
}
\arguments{
  \item{x}{
Output of the \code{sfplsim.kernel.fit} or \code{sfplsim.kNN.fit} functions (i.e. an object of the class \code{sfplsim.kernel} or \code{sfplsim.kNN}).
}
  \item{\dots}{
Further arguments.
}
 \item{object}{
Output of the \code{sfplsim.kernel.fit} or \code{sfplsim.kNN.fit} functions (i.e. an object of the class \code{sfplsim.kernel} or \code{sfplsim.kNN}).
}
}

\value{
\itemize{
\item{The matched call.}
\item{The optimal value of the tunning parameter (\code{h.opt} or \code{k.opt}).}
\item{Coefficients of \eqn{\hat{\theta}} in the B-spline basis (\code{theta.est}): a vector of \code{length(order.Bspline+nknot.theta).}}
\item{The estimated vector of linear coefficients (\code{beta.est}).}
\item{The number of non-zero components in \code{beta.est}.}
\item{The indexes of the non-zero components in \code{beta.est}.}
\item{The optimal value of the penalisation parameter (\code{lambda.opt}).}
\item{The optimal value of the criterion function, i.e. the value  obtained with \code{lambda.opt}, \code{vn.opt}  and \code{h.opt}/\code{k.opt}}
\item{Minimum value of the penalised least-squares function. That is, the value obtained using \code{theta.est}, \code{beta.est} and \code{lambda.opt}.}
\item{The penalty function used.}
\item{The criterion used to select the tuning parameter, the penalisation parameter and \code{vn}.}
\item{The optimal value of \code{vn}.}

}
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
\code{sfplsim.kernel.fit} and \code{sfplsim.kNN.fit}.
}

