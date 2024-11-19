\name{plotnls}
\alias{plotnls}

\title{Plotting diagnostics for an \code{nls} object}

\description{
This function is aimed at providing some types of plots to assess the goodness of fit for the selected model. Three plots (selectable by the argument 'which') are currently available: a plot of residuals against fitted values (which = 1), a Normal Q-Q plot (which = 2) and a plot of predicted against expected (line) and observed (symbols). By default, type = 3 is provided. As for the third graph, we can either plot all the data (type= "all") or the group means (type = "means"; the default)
}

\usage{
  plotnls(x, type = c("average", "all"),
                     xlim = NULL, gridsize = 100,
                     which = 3, ...)
}

\arguments{
  \item{x}{an object of class 'nls'}
  \item{type}{it can be either "means" or "all". In the first case, the group means are plotted for the third graph. It is only considered when which = 3}
  \item{xlim}{The limits for the x-axis (x1, x2)}
  \item{gridsize}{For 'which = 3', it sets the resolution of the fitted line}
  \item{which}{The type of graph: can be 1, 2 or 3 (see description). It defaults to 3.}
  \item{...}{additional graphical arguments}
}

\value{No return value, it produces a plot}

\details{
  It mimicks the behaviour of the function plot.lm(). It has been named 'plotnls' in order to avoid
  conflicts with the 'plot.nls' method in the 'nlme' library
}

\author{
  Andrea Onofri
}

\examples{
library(statforbiology)
degradation <- read.csv("https://www.casaonofri.it/_datasets/degradation.csv")
mod <- nls(Conc ~ A*exp(-k*Time),
               start=list(A=100, k=0.05),
               data=degradation)
plotnls(mod, which = 3)
}

