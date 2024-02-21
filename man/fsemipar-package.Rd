\name{fsemipar-package}
\alias{fsemipar-package}
\alias{fsemipar}
\docType{package}
\title{
\packageTitle{fsemipar}
}
\description{ This package is devoted to estimation or simultaneous estimation and variable selection of several functional semiparametric models with scalar response, such as the functional single-index model, the semi-functional partial linear model or the semi-functional partial linear single-index model. It also contains algorithms for addressing estimation and variable selection in the linear model, the multi-functional partial linear model, and the multi-functional partial linear single-index model when the scalar covariates with linear effects come from the discretization of a curve. In addition, the package has routines for kernel- and kNN-based estimation with Nadaraya-Watson weights of models with a nonparametric component. It also contains functions to compute predictions from all the considered models and estimation procedures.
}
\details{
The package could be divided in several tematic sections:
\enumerate{
\item Estimation and prediction of the functional single-index model.
\itemize{
\item \code{\link{projec}}.
\item \code{\link{semimetric.projec}}.
\item \code{\link{fsim.kernel.fit}} and \code{\link{fsim.kNN.fit}}.
\item \code{\link{fsim.kernel.test}} and \code{\link{fsim.kNN.test}}.
\item \code{predict, summary} and \code{print} methods for \code{fsim.kernel} and \code{fsim.kNN} classes.
}
\item Variable selection, estimation and prediction of the semi-functional partial linear single-index model.
\itemize{
\item \code{\link{sfplsim.kernel.fit}} and \code{\link{sfplsim.kNN.fit}}.
\item \code{predict, summary} and \code{print} methods for \code{sfplsim.kernel} and \code{sfplsim.kNN} classes.
}
\item Variable selection, estimation and prediction of the semi-functional partial linear model.
\itemize{
\item \code{\link{sfpl.kernel.fit}} and \code{\link{sfpl.kNN.fit}}.
\item \code{predict, summary} and \code{print} methods for \code{sfpl.kernel} and \code{sfpl.kNN} classes.
}
\item Variable selection, estimation and prediction of the linear model.
\itemize{
\item \code{\link{lm.pels.fit}}.
\item \code{predict, summary} and \code{print} methods for \code{lm.pels} class.
}

\item Variable selection, estimation and prediction of the linear model with covariates coming from the discretisation of a curve.
\itemize{
\item \code{\link{PVS.fit}}.
\item \code{predict, summary} and \code{print} methods for \code{PVS} class.
}
\item Variable selection, estimation and prediction of the multi-functional partial linear model.
\itemize{
\item \code{\link{PVS.kernel.fit}} and \code{\link{PVS.kNN.fit}}.
\item \code{predict, summary} and \code{print} methods for \code{PVS.kernel} and \code{PVS.kNN} classes.
}
\item Variable selection, estimation and prediction of the multi-functional partial linear single-index model.
\itemize{
\item \code{\link{FASSMR.kernel.fit}} and \code{\link{FASSMR.kNN.fit}}.
\item \code{\link{IASSMR.kernel.fit}} and \code{\link{IASSMR.kNN.fit}}.
\item \code{predict, summary} and \code{print} methods for \code{FASSMR.kernel}, \code{FASSMR.kNN}, \code{IASSMR.kernel} and \code{IASSMR.kNN} classes.
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


