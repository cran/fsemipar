\name{print.summary.lm}
\alias{print.lm.pels}
\alias{print.PVS}
\alias{summary.lm.pels}
\alias{summary.PVS}
\title{
Summarise information from linear models estimation
}
\description{
\code{summary} and \code{print} functions for \code{lm.pels.fit} and \code{PVS.fit}.
}
\usage{
\method{print}{lm.pels}(x, \dots)
\method{print}{PVS}(x, \dots)
\method{summary}{lm.pels}(object, \dots)
\method{summary}{PVS}(object, \dots)
}
\arguments{
  \item{x}{
Output of the \code{lm.pels.fit} or \code{PVS.fit} functions (i.e. an object of the class \code{lm.pels} or \code{PVS}).
}
  \item{\dots}{Further arguments.
}
\item{object}{
Output of the \code{lm.pels.fit} or \code{PVS.fit} functions (i.e. an object of the class \code{lm.pels} or \code{PVS}).
}
}

\value{
\itemize{
\item{The matched call.}
\item{The estimated intercept of the model.}
\item{The estimated vector of linear coefficients (\code{beta.est}).}
\item{The number of non-zero components in \code{beta.est}.}
\item{The indexes of the non-zero components in \code{beta.est}.}
\item{The optimal value of the penalisation parameter (\code{lambda.opt}).}
\item{The optimal value of the criterion function, i.e. the value  obtained with \code{lambda.opt} and \code{vn.opt} (and \code{w.opt} in the case of the PVS).}
\item{Minimum value of the penalised least-squares function. That is, the value obtained using \code{beta.est} and \code{lambda.opt}.}
\item{The penalty function used.}
\item{The criterion used to select the penalisation parameter and \code{vn}.}
\item{The optimal value of \code{vn} in the case of the \code{lm.pels} object.}
}

In the case of the \code{PVS} objects, these functions also return
the optimal number of covariates required to construct the reduced model in the first step of the algorithm (\code{w.opt}). This value is selected using the same criterion employed for selecting the penalisation parameter.

}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}

\seealso{
\code{\link{lm.pels.fit}} and \code{\link{PVS.fit}}.
}

