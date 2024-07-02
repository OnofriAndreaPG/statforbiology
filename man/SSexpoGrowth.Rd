\name{expoGrowth}
\alias{expoGrowth.fun}
\alias{NLS.expoGrowth}
\alias{DRC.expoGrowth}

\title{Exponential growth function}

\description{
  These functions provide the exponential growth equation (expoGrowth), the exponential growth equation with self-starter for the \code{\link{nls}} function and the exponential growth equation with self-starter for the \code{\link[drc]{drm}} function in the drc package.
}

\usage{
  expoGrowth.fun(predictor, init, k)
  NLS.expoGrowth(predictor, init, k)
  DRC.expoGrowth(fixed = c(NA, NA), names = c("init", "k"))
}

\arguments{
\item{predictor}{a numeric vector of values at which to evaluate the model.}
\item{init}{model parameter: response at predictor = 0}
\item{k}{model parameter: growth rate}
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The exponential growth is given by the following function:
  \deqn{ f(x) = \textrm{C0} \cdot \exp ( k x )\}}
  }

\value{
  expoGrowth.fun() and NLS.expoGrowth() return a numeric value, while DRC.expoGrowth() returns a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
}

\author{Andrea Onofri}

\note{DRC.expoGrowth() is for use with the function \code{\link[drc]{drm}}.}

\examples{
library(statforbiology)
Time <- 1:20
Resp <- c(0.18, 0.64, 1.14, 0.67, 0.32, 0.86, 0.70, 0.73, 0.89, 0.48,
    2.20, 1.03, 1.14, 2.14, 1.31, 2.08, 1.85, 1.47, 1.98, 1.30)
model <- drm(Resp ~ Time, fct = DRC.expoGrowth())
summary(model)
model2 <- nls(Resp ~ NLS.expoGrowth(Time, a, b))
summary(model2)
}

\keyword{models}
\keyword{nonlinear}
