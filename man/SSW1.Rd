% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSW1.R, R/SSW2.R
\name{SSW1}
\alias{SSW1}
\alias{W1.4.fun}
\alias{W1.3.fun}
\alias{W1.2.fun}
\alias{NLS.W1.4}
\alias{NLS.W1.3}
\alias{NLS.W1.2}
\alias{W2.4.fun}
\alias{W2.3.fun}
\alias{W2.2.fun}
\alias{NLS.W2.4}
\alias{NLS.W2.3}
\alias{NLS.W2.2}
\title{Weibull equation (Type I)}
\usage{
W1.4.fun(predictor, b, c, d, e)
W1.3.fun(predictor, b, d, e)
W1.2.fun(predictor, b, e)
NLS.W1.4(predictor, b, c, d, e)
NLS.W1.3(predictor, b, d, e)
NLS.W1.2(predictor, b, e)

W2.4.fun(predictor, b, c, d, e)
W2.3.fun(predictor, b, d, e)
W2.2.fun(predictor, b, e)
NLS.W2.4(predictor, b, c, d, e)
NLS.W2.3(predictor, b, d, e)
NLS.W2.2(predictor, b, e)
}
\arguments{
\item{predictor}{a numeric vector of values at which to evaluate the model}

\item{b}{model parameter (slope at inflection point)}

\item{c}{model parameter (lower asymptote)}

\item{d}{model parameter (higher asymptote)}

\item{e}{model parameter (abscissa at inlection point)}
}
\value{
All these functions return a numeric value

All these functions return a numeric value
}
\description{
These functions provide the Weibull equation (type I), that has an
asymmetric sygmoidal shape and it has been used for bioassay work.
These functions provide the 4-, 3- and 2-parameter equations
(W1.4.fun(), W1.3.fun() and W1.2.fun()) as well as the self-starters
for the \code{\link{nls}} function (NLS.W1.4(), NLS.W1.3() and NLS.W1.2()

These functions provide the Weibull equation (type II), that has an
asymmetric sygmoidal shape and it has been used for bioassay work.
These functions provide the 4-, 3- and 2-parameter equations
(W2.4.fun(), W2.3.fun() and W2.2.fun()) as well as the self-starters
for the \code{\link{nls}} function (NLS.W2.4(), NLS.W2.3() and NLS.W2.2()
}
\details{
These functions provide the Weibull (Type I) equation for bioassay work
This equation (4-parameters) is parameterised as:
\deqn{ f(x) = c + (d - c) \exp ( - \exp ( - b\,(\log(x) - \log(e)))) }
For the 3- and 2-parameters model, c is equal to 0, while for the 2-parameter
 model d is equal to 1.

These functions provide the Weibull (Type I) equation for bioassay work
This equation (4-parameters) is parameterised as:
\deqn{ f(x) = c + (d - c) (1 - \exp( - \exp (b \, (\log(x) - \log(e))))) }
For the 3- and 2-parameters model, c is equal to 0, while for the 2-parameter
 model d is equal to 1.
}
\examples{
library(statforbiology)
WeedDens <- c(0, 5, 10, 20, 25)
YieldLoss <- c(0, 17.9, 21.5, 27.4, 29.5)


library(statforbiology)
WeedDens <- c(0, 5, 10, 20, 25)
YieldLoss <- c(0, 17.9, 21.5, 27.4, 29.5)


}
\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}

Ritz, C., Jensen, S.M., Gerhard, D., Streibig, J.C., 2019. Dose-response analysis using R, CRC Press. ed. USA.

Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}

Ritz, C., Jensen, S.M., Gerhard, D., Streibig, J.C., 2019. Dose-response analysis using R, CRC Press. ed. USA.
}
\author{
Andrea Onofri
}