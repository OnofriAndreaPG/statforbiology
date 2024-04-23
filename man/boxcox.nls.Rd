\name{boxcox.nls}
\alias{boxcox.nls}
\alias{summary.nlsbc}

\title{Transform-both-sides (TBS) method for nonlinear regression}

\description{
  Finds the optimal Box-Cox transformation for non-linear regression models.
}

\details{
  \code{boxcox.nls} is very similar to the \code{\link[MASS]{boxcox}} in its
  arguments.
  The optimal lambda value is determined using a profile likelihood approach: 
  For each lambda value the non-linear regression model is fitted and the lambda
  value resulting in thre largest value of the log likelihood function is picked.
  If a self starter model was used in the model fit, then gradient information
  will be used in the profiling.
}

\usage{
  \method{boxcox}{nls}(object, lambda = seq(-2, 2, 1/10), plotit = TRUE, 
             start, eps = 1/50, bcAdd = 0, level = 0.95, 
             xlab = expression(lambda), ylab = "log-likelihood", ...)
  \method{summary}{nlsbc}(object, ...)
}

\arguments{
  \item{object}{object of class \code{\link{nls}}. For \code{bcSummary} the \code{nls} fit should have been obtained using \code{boxcox.nls}}
  \item{lambda}{numeric vector of lambda values; the default is (-2, 2) in steps of 0.1.}
  \item{plotit}{logical which controls whether the result should be plotted.}
  \item{start}{a list of starting values (optional).}
  \item{eps}{numeric value: the tolerance for lambda = 0; defaults to 0.02.}
  \item{bcAdd}{numeric value specifying the constant to be added on both sides prior to Box-Cox transformation. 
  The default is 0.}
  \item{level}{numeric value: the confidence level required.}
  \item{xlab}{character string: the label on the x axis, defaults to "lambda".}
  \item{ylab}{character string: the label on the y axis, defaults to "log-likelihood".}
  \item{\dots}{additional graphical parameters.}
}

\value{
  An object of class \code{\link{nls}} (returned invisibly). 
  If plotit = TRUE a plot of loglik vs lambda is shown indicating a confidence interval (by default 95\%) about 
  the optimal lambda value.
}

\references{
  Carroll, R. J. and Ruppert, D. (1988) \emph{Transformation and Weighting in Regression}, New York: Chapman and Hall (Chapter 4).
}

\author{Christian Ritz, modified by Andrea Onofri}

%\note{}

\seealso{For linear regression the analogue is \code{\link[MASS]{boxcox}}.}

\examples{
## Fitting log-logistic model without transformation
ryegrass.m1 <- nls(rootl ~ NLS.L4(conc, b, c, d, e), 
                 data=ryegrass)
summary(ryegrass.m1)

## Fitting the same model with optimal Box-Cox transformation
ryegrass.m2 <- boxcox(ryegrass.m1, plotit = TRUE)
summary(ryegrass.m2)
}

\keyword{models}
\keyword{nonlinear}

\concept{Box-Cox transform-both-sides transformation}
