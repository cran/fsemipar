\name{print.summary.sfpl}
\alias{print.sfpl.kernel}
\alias{print.sfpl.kNN}
\alias{summary.sfpl.kernel}
\alias{summary.sfpl.kNN}
\title{
Summarize information from semi-functional partial linear model (SFPLM) estimation methods
}
\description{
\code{summary} and \code{print} functions for \code{sfpl.kNN.fit} and \code{sfpl.kernel.fit}.
}
\usage{
\method{print}{sfpl.kernel}(x, \dots)
\method{print}{sfpl.kNN}(x, \dots)
\method{summary}{sfpl.kernel}(object, \dots)
\method{summary}{sfpl.kNN}(object, \dots)
}
\arguments{
  \item{x}{
Output of the \code{sfpl.kernel.fit} or \code{sfpl.kNN.fit} functions (i.e. an object of the class \code{sfpl.kernel} or \code{sfpl.kNN}).
}
  \item{\dots}{
Further arguments.
}
\item{object}{
Output of the \code{sfpl.kernel.fit} or \code{sfpl.kNN.fit} functions (i.e. an object of the class \code{sfpl.kernel} or \code{sfpl.kNN}).
}
}

\value{
\itemize{
\item{The matched call.}
\item{The optimal value of the tunning parameter (\code{h.opt} or \code{k.opt}).}
\item{The estimated vector of linear coefficients (\code{beta.est}).}
\item{The number of non-zero components in \code{beta.est}.}
\item{The indexes of the non-zero components in \code{beta.est}.}
\item{The optimal value of the penalisation parameter (\code{lambda.opt}).}
\item{The optimal value of the criterion function, i. e. the value  obtained with \code{lambda.opt}, \code{vn.opt}  and \code{h.opt}/\code{k.opt}}
\item{Minimum value of the penalized profile least-squares function. That is, the value obtained using \code{beta.est}.}
\item{The penalty function used.}
\item{The criterion used to select the tunning parameter, the penalisation parameter and \code{vn}.}
\item{The optimal value of \code{vn}.}
}
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
\code{sfpl.kernel.fit} and \code{sfpl.kNN.fit}.
}
