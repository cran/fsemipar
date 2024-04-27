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
Graphical representation of regression model outputs
}
\description{
\code{plot} functions to generate visual representations for the outputs of several fitting functions:  \code{FASSMR.kernel.fit}, \code{FASSMR.kNN.fit},  \code{fsim.kernel.fit},  \code{fsim.kernel.fit.optim}, \code{fsim.kNN.fit}, \code{fsim.kNN.fit.optim}, \code{IASSMR.kernel.fit}, \code{IASSMR.kNN.fit}, \code{lm.pels.fit}, \code{PVS.fit}, \code{PVS.kernel.fit}, \code{PVS.kNN.fit}, \code{sfpl.kernel.fit}, \code{sfpl.kNN.fit},\code{sfplsim.kernel.fit} and \code{sfplsim.kNN.fit}.
}
\usage{
\method{plot}{FASSMR.kernel}(x,ind=1:10, size=15,col1=1,col2=2,col3=4,option=0,\dots)

\method{plot}{FASSMR.kNN}(x,ind=1:10, size=15,col1=1,col2=2,col3=4,option=0, \dots)

\method{plot}{fsim.kernel}(x,size=15,col1=1,col2=2, \dots)

\method{plot}{fsim.kNN}(x,size=15,col1=1,col2=2,\dots)

\method{plot}{IASSMR.kernel}(x,ind=1:10, size=15,col1=1,col2=2,col3=4,option=0, \dots)

\method{plot}{IASSMR.kNN}(x,ind=1:10, size=15,col1=1,col2=2,col3=4,option=0, \dots)

\method{plot}{lm.pels}(x,size=15,col1=1,col2=2,col3=4, \dots)

\method{plot}{PVS}(x,ind=1:10, size=15,col1=1,col2=2,col3=4,option=0, \dots)

\method{plot}{PVS.kernel}(x,ind=1:10, size=15,col1=1,col2=2,col3=4,option=0, \dots)

\method{plot}{PVS.kNN}(x,ind=1:10, size=15,col1=1,col2=2,col3=4,option=0, \dots)

\method{plot}{sfpl.kernel}(x,size=15,col1=1,col2=2,col3=4, \dots)

\method{plot}{sfpl.kNN}(x,size=15,col1=1,col2=2,col3=4, \dots)

\method{plot}{sfplsim.kernel}(x,size=15,col1=1,col2=2,col3=4, \dots)

\method{plot}{sfplsim.kNN}(x,size=15,col1=1,col2=2,col3=4, \dots)
}

\arguments{
  \item{x}{
Output of the functions mentioned in the \code{Description} (i.e. an object of the class \code{FASSMR.kernel}, \code{FASSMR.kNN}, \code{fsim.kernel},\code{fsim.kNN}, \code{IASSMR.kernel}, \code{IASSMR.kNN}, \code{lm.pels}, \code{PVS}, \code{PVS.kernel}, \code{PVS.kNN}, \code{sfpl.kernel},\code{sfpl.kNN}, \code{sfplsim.kernel} or \code{sfplsim.kNN}).
}
\item{ind}{Indexes of the colors for the curves in the chart of estimated impact points. The default is \code{1:10}}
\item{size}{
The size for title and axis labels in pts. The default is 15.
}
  \item{col1}{
Color of the points in the charts. Also, color of the estimated functional index representation. The default is black.
}
  \item{col2}{
Color of the nonparametric fit representation in FSIM functions, and of the straight line in 'Response vs Fitted Values' charts. The default is red.
}
  \item{col3}{
Color of the nonparametric fit of the residuals in 'Residuals vs Fitted Values' charts. The default is blue.
}
\item{option}{Selection of charts to be plotted. The default, \code{option = 0}, means all charts are plotted. See the section \code{Details}.}
\item{\dots}{Further arguments passed to or from other methods.

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
\item For the classes \code{lm.pels}, \code{sfpl.kernel} and \code{sfpl.kNN}:
\enumerate{
\item The response over the \code{fitted.values}.
\item The \code{residuals} over the \code{fitted.values}.
}

\item For the classes \code{sfplsim.kernel} and \code{sfplsim.kNN}:
\enumerate{
     \item  The estimated functional index: \eqn{\hat{\theta}}.
     \item The response over the \code{fitted.values}.
     \item The \code{residuals} over the \code{fitted.values}.
 }
\item For the classes \code{FASSMR.kernel}, \code{FASSMR.kNN}, \code{IASSMR.kernel}, \code{IASSMR.kNN}, \code{sfplsim.kernel} and \code{sfplsim.kNN}:
\enumerate{
    \item If \code{option=1}: The curves with the estimated impact points (in dashed vertical lines).
     \item If \code{option=2}: The estimated functional index: \eqn{\hat{\theta}}.
     \item If \code{option=3}:
     \itemize{
     \item The response over the \code{fitted.values}.
     \item The \code{residuals} over the \code{fitted.values}.}
     \item If \code{option=0}: All chart are plotted.
 }

\item For the classes \code{PVS}, \code{PVS.kNN}, \code{PVS.kernel}:
\enumerate{
    \item If \code{option=1}: The curves with the estimated impact points (in dashed vertical lines).
     \item If \code{option=2}: \itemize{ 
     \item The response over the \code{fitted.values}.
     \item The \code{residuals} over the \code{fitted.values}.}
     \item If \code{option=0}: All chart are plotted.
 }
}
All the routines implementing the plot S3 method use internally the R package \code{ggplot2} to produce elegant and high quality
charts.

}

\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}



\seealso{
\code{\link{FASSMR.kernel.fit}}, \code{\link{FASSMR.kNN.fit}},  \code{\link{fsim.kernel.fit}}, \code{\link{fsim.kNN.fit}}, \code{\link{IASSMR.kernel.fit}}, \code{\link{IASSMR.kNN.fit}}, \code{\link{lm.pels.fit}}, \code{\link{PVS.fit}}, \code{\link{PVS.kernel.fit}}, \code{\link{PVS.kNN.fit}}, \code{\link{sfpl.kernel.fit}}, \code{\link{sfpl.kNN.fit}}, \code{\link{sfplsim.kernel.fit}} and \code{\link{sfplsim.kNN.fit}}.
}

