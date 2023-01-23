\name{Sugar}
\alias{Sugar}
\docType{data}
\title{
Sugar data
}
\description{
Ash content and absorbance spectra at two different excitation wavelengths of 268 samples of sugar. Detailed information about this dataset can be found in \url{https://ucphchemometrics.com/datasets/}.
}
\usage{data(Sugar)}
\format{
  A \code{list} containing:
 \itemize{
 \item{\code{ash}: A vector with the ash contents.}
 \item{\code{wave.290}: A matrix containing the absorbance spectra observed on 571 equally spaced wavelengths in the range 275-560 nm at excitation wavelengths 290 nm.}
 \item{\code{wave.240}: A matrix containing the absorbance spectra observed on 571 equally spaced wavelengths in the range 275-560 nm at excitation wavelengths 240 nm.}
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
