\name{fsemipar-package}
\alias{fsemipar-package}
\alias{fsemipar}
\docType{package}
\title{
\packageTitle{fsemipar}
}
\description{This package is dedicated to the estimation and simultaneous estimation and variable selection in several functional semiparametric models with scalar response. These include the functional single-index model, the semi-functional partial linear model, and the semi-functional partial linear single-index model. Additionally, it encompasses algorithms for addressing estimation and variable selection in linear models and bi-functional partial linear models  when the scalar covariates with linear effects are derived from the discretisation of a curve. Furthermore, the package offers routines for kernel- and kNN-based estimation using Nadaraya-Watson weights in models with a nonparametric or semiparametric component. It also includes S3 methods (predict, plot, print, summary) to facilitate statistical analysis across all the considered models and estimation procedures.
}
\details{
The package can be divided into several thematic sections:
\enumerate{
\item Estimation of the functional single-index model.
\itemize{
\item \code{\link{projec}}.
\item \code{\link{semimetric.projec}}.
\item \code{\link{fsim.kernel.fit}} and \code{\link{fsim.kNN.fit}}.
\item \code{\link{fsim.kernel.fit.optim}} and \code{\link{fsim.kNN.fit.optim}}
\item \code{\link{fsim.kernel.test}} and \code{\link{fsim.kNN.test}}.
\item \code{predict, plot, summary} and \code{print} methods for \code{fsim.kernel} and \code{fsim.kNN} classes.
}
\item Simultaneous estimation and variable selection in linear and semi-functional partial linear models.
\enumerate{
\item Linear model
\itemize{
\item \code{\link{lm.pels.fit}}. 
\item \code{predict, summary, plot} and \code{print} methods for \code{lm.pels} class.
}
\item Semi-functional partial linear model.
\itemize{
\item \code{\link{sfpl.kernel.fit}} and \code{\link{sfpl.kNN.fit}}.
\item \code{predict, summary, plot} and \code{print} methods for \code{sfpl.kernel} and \code{sfpl.kNN} classes.
}
\item Semi-functional partial linear single-index model.
\itemize{
\item \code{\link{sfplsim.kernel.fit}} and \code{\link{sfplsim.kNN.fit}}.
\item \code{predict, summary, plot} and \code{print} methods for \code{sfplsim.kernel} and \code{sfplsim.kNN} classes.
}
}
\item Algorithms for impact point selection in models with covariates derived from the discretisation of a curve.
\enumerate{
\item Linear model
\itemize{
\item \code{\link{PVS.fit}}. 
\item \code{predict, summary, plot} and \code{print} methods for \code{PVS} class.
}
\item Bi-functional partial linear model.
\itemize{
\item \code{\link{PVS.kernel.fit}} and \code{\link{PVS.kNN.fit}}.
\item \code{predict, summary, plot} and \code{print} methods for \code{PVS.kernel} and \code{PVS.kNN} classes.
}
\item Bi-functional partial linear single-index model.
\itemize{
\item \code{\link{FASSMR.kernel.fit}} and \code{\link{FASSMR.kNN.fit}}.
\item \code{\link{IASSMR.kernel.fit}} and \code{\link{IASSMR.kNN.fit}}.
\item \code{predict, summary, plot} and \code{print} methods for \code{FASSMR.kernel}, \code{FASSMR.kNN}, \code{IASSMR.kernel} and \code{IASSMR.kNN} classes.
}
}
\item Two datasets: \code{\link{Tecator}} and \code{\link{Sugar}}.
}


}
\author{
\packageAuthor{fsemipar}

Maintainer: \packageMaintainer{fsemipar}
}
\references{

Aneiros, G. and Vieu, P., (2014) Variable selection in infinite-dimensional problems, \emph{Statistics and Probability Letters}, \bold{94}, 12--20. \doi{https://doi.org/10.1016/j.spl.2014.06.025}.

Aneiros, G., Ferraty, F., and Vieu, P., (2015) Variable selection in partial linear regression with functional
covariate, \emph{Statistics}, \bold{49} 1322--1347, \doi{https://doi.org/10.1080/02331888.2014.998675}.

Aneiros, G., and Vieu, P., (2015) Partial linear modelling with multi-functional covariates. \emph{Computational Statistics}, \bold{30}, 647--671. \doi{https://doi.org/10.1007/s00180-015-0568-8}.

Novo S., Aneiros, G., and Vieu, P., (2019) Automatic and location-adaptive estimation in functional single-index regression, \emph{Journal of Nonparametric Statistics}, \bold{31(2)}, 364--392, \doi{https://doi.org/10.1080/10485252.2019.1567726}.

Novo, S., Aneiros, G., and Vieu, P., (2021) Sparse semiparametric regression
when predictors are mixture of functional and high-dimensional variables, \emph{TEST},
\bold{30}, 481--504, \doi{https://doi.org/10.1007/s11749-020-00728-w}.

Novo, S., Aneiros, G., and Vieu, P., (2021) A kNN procedure in semiparametric
functional data analysis, \emph{Statistics and Probability Letters}, \bold{171}, 109028, \doi{https://doi.org/10.1016/j.spl.2020.109028}.

Novo, S., Vieu, P., and Aneiros, G., (2021) Fast and efficient algorithms for
sparse semiparametric bi-functional regression, \emph{Australian and New Zealand
Journal of Statistics}, \bold{63}, 606--638, \doi{https://doi.org/10.1111/anzs.12355}.
}


