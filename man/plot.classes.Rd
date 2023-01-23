\name{plot.classes}

\alias{plot.FASSMR.kernel}
\alias{plot.FASSMR.kNN}
\alias{plot.fsim.kernel}
\alias{plot.fsim.kNN}
\alias{plot.IASSMR.kernel}
\alias{plot.IASSMR.kNN}
\alias{plot.lm.pels}
\alias{plot.PVS}
\alias{plot.PVS.kernel}
\alias{plot.PVS.kNN}
\alias{plot.sfpl.kernel}
\alias{plot.sfpl.kNN}
\alias{plot.sfplsim.kernel}
\alias{plot.sfplsim.kNN}


\title{
Plot outputs from regression estimation methods
}
\description{
\code{plot} function for \code{FASSMR.kernel.fit}, \code{FASSMR.kNN.fit},  \code{fsim.kernel.fit}, \code{fsim.kNN.fit}, \code{IASSMR.kernel.fit}, \code{IASSMR.kNN.fit}, \code{lm.pels.fit}, \code{PVS.fit}, \code{PVS.kernel.fit}, \code{PVS.kNN.fit}, \code{sfpl.kernel.fit}, \code{sfpl.kNN.fit},\code{sfplsim.kernel.fit} and \code{sfplsim.kNN.fit}.
}
\usage{
\method{plot}{FASSMR.kernel}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{FASSMR.kNN}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{fsim.kernel}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{fsim.kNN}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{IASSMR.kernel}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{IASSMR.kNN}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{lm.pels}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{PVS}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{PVS.kernel}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{PVS.kNN}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{sfpl.kernel}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{sfpl.kNN}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{sfplsim.kernel}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)

\method{plot}{sfplsim.kNN}(x, cex.axis = 1.5, cex.lab = 1.5, cex = 2, col = 1, cex.main = 1.5, \dots)
}

\arguments{
  \item{x}{
Output of the functions mentioned in the \code{Description} (i.e. an object of the class \code{FASSMR.kernel}, \code{FASSMR.kNN}, \code{fsim.kernel},\code{fsim.kNN}, \code{IASSMR.kernel}, \code{IASSMR.kNN}, \code{lm.pels}, \code{PVS}, \code{PVS.kernel}, \code{PVS.kNN}, \code{sfpl.kernel},\code{sfpl.kNN}, \code{sfplsim.kernel} or \code{sfplsim.kNN}).
}
  \item{cex.axis}{
The magnification to be used for axis annotation relative to the current setting of \code{cex}. The default is 1.5.
}
  \item{cex.lab}{
The magnification to be used for x and y labels relative to the current setting of \code{cex}. The default is 1.5.
}
  \item{cex}{
A numerical value giving the amount by which plotting text and symbols should be magnified. The default is 2.
}
  \item{col}{
A specification for the default plotting color. The default is \code{color=1}.
}
  \item{cex.main}{
The magnification to be used for main titles relative to the current setting of \code{cex}. The default is 1.5.
}
  \item{\dots}{

}
}

\value{
The functions return different graphical representations.
\itemize{
\item For the classes \code{fsim.kNN} and \code{fsim.kernel}:
     \enumerate{
     \item The estimated functional index: \eqn{\hat{\theta}}.
     \item The regression fit.
     }

\item For the classes \code{FASSMR.kernel}, \code{FASSMR.kNN}, \code{IASSMR.kernel}, \code{IASSMR.kNN}, \code{sfplsim.kernel} and \code{sfplsim.kNN}.
\enumerate{
     \item The response over the \code{fitted.values}.
      \item The \code{residuals} over the \code{fitted.values}.
     \item The estimated functional index: \eqn{\hat{\theta}}.
     }

\item For the classes \code{lm.pels}, \code{PVS}, \code{PVS.kernel}, \code{PVS.kNN}, \code{sfpl.kernel} and \code{sfpl.kNN}.
\enumerate{
\item The response over the \code{fitted.values}.
\item The \code{residuals} over the \code{fitted.values}.
}
}
}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}



\seealso{
\code{\link{FASSMR.kernel.fit}}, \code{\link{FASSMR.kNN.fit}},  \code{\link{fsim.kernel.fit}}, \code{\link{fsim.kNN.fit}}, \code{\link{IASSMR.kernel.fit}}, \code{\link{IASSMR.kNN.fit}}, \code{\link{lm.pels.fit}}, \code{\link{PVS.fit}}, \code{\link{PVS.kernel.fit}}, \code{\link{PVS.kNN.fit}}, \code{\link{sfpl.kernel.fit}}, \code{\link{sfpl.kNN.fit}}, \code{\link{sfplsim.kernel.fit}} and \code{\link{sfplsim.kNN.fit}}.
}

