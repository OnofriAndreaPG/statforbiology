\name{expoDecay}
\alias{expoDecay.fun}
\alias{NLS.expoDecay}
\alias{DRC.expoDecay}

\title{Exponential decay function}

\description{
  These functions provide the exponential decay function (expoDecay), the exponential decay function with self-starter for the \code{\link{nls}} function and the exponential decay function with self-starter for the \code{\link[drc]{drm}} function in the drc package.
}

\usage{
  expoDecay.fun(predictor, C0, k)
  NLS.expoDecay(predictor, C0, k)
  DRC.expoDecay(fixed = c(NA, NA), names = c("C0", "k"))
}

\arguments{
\item{predictor}{a numeric vector of values at which to evaluate the model.}
\item{C0}{model parameter}
\item{k}{model parameter}
 \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
 \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
}

\details{
  The exponential decay is given by the following function:
  \deqn{ f(x) = \textrm{C0} \cdot \exp (- k x )\}}
  }

\value{
  expoDecay.fun and NLS.expoDecay return a numeric value, while DRC.expoDecay returns a list containing the nonlinear function, the self starter function and the parameter names.
}

\references{
Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.

Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
}

\author{Andrea Onofri}

\note{DRC.expoDecay is for use with the function \code{\link[drc]{drm}}.}

\examples{
data("degradation")
model <- drm(Conc ~ Time, fct = DRC.expoDecay(),
             data = degradation)
summary(model)
model2 <- nls(Conc ~ NLS.expoDecay(Time, a, b),
              data = degradation)
summary(model2)
}

\keyword{models}
\keyword{nonlinear}
