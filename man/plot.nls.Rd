\name{plot.nls}
\alias{plot.nls}

\title{Plotting diagnostics for an \code{nls} object}

\description{
This function is aimed at providing some types of plots to assess the goodness of fit for the selected model. Three plots (selectable by the argument 'which') are currently available: a plot of residuals against fitted values (which = 1), a Normal Q-Q plot (which = 2) and a plot of predicted against expected (line) and observed (symbols). By default, type = 3 is provided. As for the third graph, we can either plot all the data (type= "all") or the group means (type = "means"; the default)
}

\usage{

  \method{plot}{nls}(x, type = c("average", "all"),
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

\details{
  It mimicks the behaviour of the function plot.lm()
}

\author{
  Andrea Onofri
}

\examples{
degradation <- read.csv("https://www.casaonofri.it/_datasets/degradation.csv")
mod <- nls(Conc ~ A*exp(-k*Time),
               start=list(A=100, k=0.05),
               data=degradation)
plot(mod, which = 3)
}

