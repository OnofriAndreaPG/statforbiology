\name{plotRes}
\alias{plotRes}

\title{Graphical analyses of residuals for \code{drc} objects}

\description{
Two graphs of model residualsts (selectable by the argument 'which') are currently available: a plot of residuals against fitted values (which = 1) and a Normal Q-Q plot (which = 2). By default, the first graph is provided.
}

\usage{
 plotRes(x, which = 1, ...)
}

\arguments{
  \item{x}{an object of class 'drc'}
  \item{which}{The type of graph: can be 1 or 2}
  \item{...}{additional graphical arguments}
}

\details{
  It mimicks the behaviour of the function plot.lm()
}

\author{
  Andrea Onofri
}

\examples{
library(statforbiology)
data(degradation)
mod <- nls(Conc ~ NLS.expoDecay(Time, a, k),
               data = degradation)
plot(mod, which = 2)
}

