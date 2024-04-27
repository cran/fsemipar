\name{semimetric.projec}
\alias{semimetric.projec}
\title{
Projection semi-metric computation
}
\description{
Computes the projection semi-metric between each curve in \code{data1} and each curve in \code{data2}, given a functional index \eqn{\theta}.
}
\usage{
semimetric.projec(data1, data2, theta, order.Bspline = 3, nknot.theta = 3,
  range.grid = NULL, nknot = NULL)
}
\arguments{
  \item{data1}{
Matrix containing functional data collected by row.
}
  \item{data2}{
Matrix containing functional data collected by row.
}
  \item{theta}{
Vector containing the coefficients of \eqn{\theta} in a B-spline basis, so that \code{length(theta)=order.Bspline+nknot.theta}.
}
  \item{order.Bspline}{
Order of the B-spline basis functions for the B-spline representation of \eqn{\theta}. This is the number of coefficients in each piecewise polynomial segment. The default is 3.
}
  \item{nknot.theta}{
Number of regularly spaced interior knots of the B-spline basis. The default is 3.
}
  \item{range.grid}{
Vector of length 2 containing  the range of the discretisation of the functional data. If \code{range.grid=NULL}, then \code{range.grid=c(1,p)} is considered, where \code{p} is the discretization size of \code{data} (i.e. \code{ncol(data)}).
}
  \item{nknot}{
Number of  regularly spaced interior knots for the B-spline representation of the functional data. The default value is \code{(p - order.Bspline - 1)\%/\%2}.
}
}
\details{
For \eqn{x_1,x_2 \in \mathcal{H}, }, where \eqn{\mathcal{H}} is a separable Hilbert space, the projection semi-metric in the direction \eqn{\theta\in \mathcal{H}} is defined as   \deqn{d_{\theta}(x_1,x_2)=|\langle\theta,x_1-x_2\rangle|.}

The function \code{semimetric.projec} computes this projection semi-metric using the B-spline representation of the curves and \eqn{\theta}. The dimension of the B-spline basis for \eqn{\theta} is determined by \code{order.Bspline}+\code{nknot.theta}.

}
\value{
A matrix containing the projection semi-metrics for each pair of curves.
}
\references{
Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single--index regression. \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.
}
\author{
German Aneiros Perez \email{german.aneiros@udc.es} 

Silvia Novo Diaz  \email{snovo@est-econ.uc3m.es}
}



\seealso{
See also \code{\link{projec}}.
}
\examples{

data("Tecator")
names(Tecator)
y<-Tecator$fat
X<-Tecator$absor.spectra

#length(theta)=6=order.Bspline+nknot.theta 
semimetric.projec(data1=X[1:5,], data2=X[5:10,],theta=c(1,0,0,1,1,-1),
  nknot.theta=3,nknot=20,range.grid=c(850,1050))

}

