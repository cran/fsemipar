\name{print.summary.mfpl}
\alias{print.PVS.kernel}
\alias{print.PVS.kNN}
\alias{summary.PVS.kernel}
\alias{summary.PVS.kNN}
\title{
Summarize information from multi-functional partial linear  model (MFPLM) estimation methods
}
\description{
\code{summary} and \code{print} functions for \code{PVS.kernel.fit} and \code{PVS.kNN.fit}.
}
\usage{
\method{print}{PVS.kernel}(x, \dots)
\method{print}{PVS.kNN}(x, \dots)
\method{summary}{PVS.kernel}(object, \dots)
\method{summary}{PVS.kNN}(object, \dots)
}
\arguments{
\item{x}{
Output of the \code{PVS.kernel.fit} or \code{PVS.kNN.fit} functions (i.e. an object of the class \code{PVS.kernel} or \code{PVS.kNN}).
}
\item{\dots}{Further arguments.
}
\item{object}{
Output of the \code{PVS.kernel.fit} or \code{PVS.kNN.fit} functions (i.e. an object of the class \code{PVS.kernel} or \code{PVS.kNN}).
}
}
\value{
\itemize{
\item{The matched call.}
\item{The optimal value of the tunning parameter (\code{h.opt} or \code{k.opt}).}
\item{The optimal initial number of covariates to build the reduced model (\code{w.opt}).}
\item{The estimated vector of linear coefficients (\code{beta.est}).}
\item{The number of non-zero components in \code{beta.est}.}
\item{The indexes of the non-zero components in \code{beta.est}.}
\item{The optimal value of the penalisation parameter (\code{lambda.opt}).}
\item{The optimal value of the criterion function, i. e. the value  obtained with \code{w.opt}, \code{lambda.opt}, \code{vn.opt}  and \code{h.opt}/\code{k.opt}}
\item{Minimum value of the penalized profile least-squares function. That is, the value obtained using \code{beta.est}.}
\item{The penalty function used.}
\item{The criterion used to select the initial number of covariates used to build the reduced model, the tunning parameter, the penalisation parameter and \code{vn}.}
}

}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}


\seealso{
\code{PVS.kernel.fit} and \code{PVS.kNN.fit}.
}


