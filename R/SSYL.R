#' Yield loss equation (Rectangular hyperbola)
#'
#' These functions provide the yield loss equation, based on
#' a rectangular hyperbola, supporting a higher asymptote and
#' no inflection points. These functions provide the
#' equation (YL.fun), the equation
#' with self-starters for the \code{\link{nls}}
#' function (NLS.YL) and equation with self-starters for
#' the \code{\link[drc]{drm}} function in the drc package (DRC.YL)
#'
#'
#' @name SSYL
#' @aliases YL.fun
#' @aliases NLS.YL
#' @aliases DRC.YL
#'
#' @usage YL.fun(predictor, i, A)
#' NLS.YL(predictor, i, A)
#' DRC.YL(fixed = c(NA, NA), names = c("i", "A"))
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param i model parameter (initial slope)
#' @param A model parameter (maximum percentage yield loss)
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.
#' @param names a vector of character strings giving the names of the parameters. The default is usually reasonable.
#'
#' @details
#' The Yield-loss equation is parameterised as:
#' \deqn{ f(x) = \frac{i \, x}{1 + (i \, x)/A}   },
#' it is convex and asymptotically increasing, while the predictor increases.
#' The response is zero when the predictor is also zero and it was mainly used
#' to describe yield losses (in percentage) due to weed competition, expressed
#' as plant density (Cousens, 1985)
#'
#' @return  YL.fun and NLS.YL return a numeric value, while DRC.YL returns
#' a list containing the nonlinear function and the self starter function
#'
#' @author Andrea Onofri
#'
#' @references Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
#' @references Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
#' @references Cousens, R., 1985. A simple model relating yield loss to weed density. Annals of Applied Biology 107, 239–252. https://doi.org/10.1111/j.1744-7348.1985.tb01567.x
#'
#' @examples
#' library(statforbiology)
#' WeedDens <- c(0, 5, 10, 20, 25)
#' YieldLoss <- c(0, 17.9, 21.5, 27.4, 29.5)
#'
#' # nls fit
#' mod.nls <- nls(YieldLoss ~ NLS.YL(WeedDens, i, A) )
#' summary(mod.nls)
#' # drm fit
#' mod.drc <- drm(YieldLoss ~ WeedDens, fct = DRC.YL() )
#' summary(mod.drc)
#'
# Yield-Loss function (rectangular hyperbola) ######################
YL.fun <- function(predictor, i, A) {
  i * predictor/(1 + i/A * predictor)
}

YL.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <- xy[, "x"]; y <- xy[, "y"]
  pseudoX <- 1 / x[x > 0]; pseudoY <- 1 / y[x > 0]
  lmFit <- lm(pseudoY ~ pseudoX)
  coefs <- coef(lmFit)
  A <- 1 / coefs[1]
  i <- 1 / coefs[2]
  value <- c(i, A)
  names(value) <- mCall[c("i", "A")]
  value
}
NLS.YL <- selfStart(YL.fun, YL.Init, parameters=c("i", "A"))

"DRC.YL" <- function(fixed = c(NA, NA), names = c("i", "A")) {
  ## Checking arguments
  numParm <- 2
  if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

  ## Fixing parameters (using argument 'fixed')
  notFixed <- is.na(fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]

  ## Defining the non-linear function
  fct <- function(x, parm) {
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm
    i <- parmMat[, 1]; A <- parmMat[, 2]
    YL.fun(x, i, A)
  }

  ## Defining self starter function
  ssfct <- function(dataf) {
    x <- dataf[, 1]
    y <- dataf[, 2]

    #regression on pseudo y values
    pseudoY <- 1 /  y[x > 0]
    pseudoX <- 1 / x [x > 0]
    coefs <- coef( lm(pseudoY ~ pseudoX) )
    A <- 1 / coefs[1]; i <- 1 / coefs[2]

    return(c(i, A)[notFixed])
  }

  ## Defining names
  pnames <- names[notFixed]

  ## Defining derivatives

  ## Defining the ED function

  ## Defining the inverse function

  ## Defining descriptive text
  text <- "Yield-Loss function (Cousens, 1985)"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

  class(returnList) <- "drcMean"
  invisible(returnList)
}
