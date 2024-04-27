\name{Sugar}
\alias{Sugar}
\docType{data}
\title{
Sugar data
}
\description{
Ash content and absorbance spectra at two different excitation wavelengths of 268 sugar samples. Detailed information about this dataset can be found at \url{https://ucphchemometrics.com/datasets/}.
}
\usage{data(Sugar)}
\format{
  A \code{list} containing:
 \itemize{
 \item{\code{ash}: A vector with the ash content.}
 \item{\code{wave.290}: A matrix containing the absorbance spectra observed at 571 equally spaced wavelengths in the range of 275-560nm, at an excitation wavelengths of 290nm.}
 \item{\code{wave.240}: A matrix containing the absorbance spectra observed at 571 equally spaced wavelengths in the range of 275-560nm, at an excitation wavelengths of 240nm.}
 }
}

\references{
Aneiros, G., and Vieu, P. (2015) Partial linear modelling with multi-functional covariates. \emph{Computational Statistics}, \bold{30}, 647--671, \doi{https://doi.org/10.1007/s00180-015-0568-8}.

Novo, S., Vieu, P., and Aneiros, G., (2021) Fast and efficient algorithms for
sparse semiparametric bi-functional regression. \emph{Australian and New Zealand
Journal of Statistics}, \bold{63}, 606--638, \doi{https://doi.org/10.1111/anzs.12355}.
}
\examples{
data(Sugar)
names(Sugar)
Sugar$ash
dim(Sugar$wave.290)
dim(Sugar$wave.240)
}
\keyword{datasets}
