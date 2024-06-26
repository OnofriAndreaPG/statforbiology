% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSE.R
\name{SSE}
\alias{SSE}
\alias{E4.fun}
\alias{E3.fun}
\alias{E2.fun}
\alias{NLS.E4}
\alias{NLS.E3}
\alias{NLS.E2}
\alias{DRC.E4}
\alias{DRC.E3}
\alias{DRC.E2}
\title{Modified Gompertz equations}
\usage{
E4.fun(predictor, b, c, d, e)
E3.fun(predictor, b, d, e)
E2.fun(predictor, b, e)
NLS.E4(predictor, b, c, d, e)
NLS.E3(predictor, b, d, e)
NLS.E2(predictor, b, e)
DRC.E4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
DRC.E3(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
DRC.E2(fixed = c(NA, NA), names = c("b", "e"))
}
\arguments{
\item{predictor}{a numeric vector of values at which to evaluate the model}

\item{b}{model parameter (slope at inflection point)}

\item{c}{model parameter (lower asymptote)}

\item{d}{model parameter (higher asymptote)}

\item{e}{model parameter (abscissa at inflection point)}

\item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}

\item{names}{names. A vector of character strings giving the names of the parameters. The default is reasonable.}
}
\value{
E4.fun, E3.fun, E2.fun, NLS.E4, NLS.E3 and NLS.E2 return a numeric value,
while DRC.E4, DRC.E3 and DRC.E2 return a list containing the nonlinear function,
the self starter function and the parameter names.
}
\description{
These functions provide the modified Gompertz equations with 4 (E4.fun), 3 (E3.fun)
and 2 (E2.fun) parameters with self-starter for the \code{\link{nls}}
function (NLS.E4, NLS.E3 and NLS.E2) and for the \code{\link[drc]{drm}} function
in the 'drc' package (DRC.E4, DRC.E3 and DRC.E2).
}
\details{
The modified Gompertz equation is parameterised as:
\deqn{ f(x) = c + (d - c) \, (1 - \exp \left[-exp( b (x - e))) \right] }
It is a sygmoidally shaped curve and it is asymmetric about its inflection
point, but the type of asymmetry is different from the Gompertz equation.
For the 3- and 2-parameter model c is equal to 0, while for the
2-parameter model d is equal to 1.
}
\examples{
data(beetGrowth)
mod3 <- nls(weightInf ~ NLS.E3(DAE, b, c, d), data = beetGrowth)
summary(mod3)
plot(mod3)

}
\author{
Andrea Onofri
}
