\name{Tecator}
\alias{Tecator}
\docType{data}
\title{
Tecator data
}
\description{
Fat, protein, and moisture content, along with absorbance spectra (including the first and second derivatives), of 215 meat samples. 
A detailed description of the data can be found at \url{http://lib.stat.cmu.edu/datasets/tecator}.
}
\usage{data(Tecator)}
\format{
 A \code{list} containing:
 \itemize{
 \item{\code{fat}: A vector with the fat content.}
 \item{\code{protein}: A vector with the protein content.}
 \item{\code{moisture}: A vector with the moisture content.}
 \item{\code{absor.spectra}: A matrix containing the near-infrared absorbance spectra observed at 100 equally spaced wavelengths in the range of 850-1050nm.}
 \item{\code{absor.spectra1}: Fist derivative of the absorbance spectra (computed using B-spline representation of the curves).}
 \item{\code{absor.spectra2}: Second derivative of the absorbance spectra (computed using B-spline representation of the curves).}
 }
 
}

\references{
Ferraty, F. and Vieu, P. (2006) \emph{Nonparametric functional data analysis}, Springer Series in Statistics, New York.
}

\examples{
data(Tecator)
names(Tecator)
Tecator$fat
Tecator$protein
Tecator$moisture
dim(Tecator$absor.spectra)
}
\keyword{datasets}
