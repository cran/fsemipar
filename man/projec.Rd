\name{projec}
\alias{projec}
\title{Inner product computation}
\description{
Computes the inner product between each curve collected in \code{data} and a particular curve  \eqn{\theta}.
}
\usage{
projec(data, theta, order.Bspline = 3, nknot.theta = 3, range.grid = NULL, nknot = NULL)
}
\arguments{
  \item{data}{
Matrix containing functional data collected by row
}
  \item{theta}{
Vector containing the coefficients of \eqn{\theta} in a B-spline basis, so that \code{length(theta)=order.Bspline+nknot.theta}
}
  \item{order.Bspline}{
Order of the B-spline basis functions for the B-spline representation of \eqn{\theta}. This is the number of coefficients in each piecewise polynomial segment. The default is 3.
}
  \item{nknot.theta}{
Number of regularly spaced interior knots of the B-spline basis. The default is 3.
}
  \item{range.grid}{
Vector of length 2 containing  the range of the discretisation of the functional data. If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the discretisation size of \code{data} (i.e. \code{ncol(data)}).
}
  \item{nknot}{
Number of regularly spaced interior knots for the B-spline representation of the functional data. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
}

\value{
A matrix containing the inner products.
}
\references{
Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single--index regression. \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}
\note{
The construction of this code is based on that by Frederic Ferraty, which is available on his website \url{https://www.math.univ-toulouse.fr/~ferraty/SOFTWARES/NPFDA/index.html}.}


\seealso{
See also \code{\link{semimetric.projec}}.
}
\examples{

data("Tecator")
names(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra

#length(theta)=6=order.Bspline+nknot.theta 
projec(X,theta=c(1,0,0,1,1,-1),nknot.theta=3,nknot=20,range.grid=c(850,1050))
}

