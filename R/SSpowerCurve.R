#' Power curve equation
#'
#' These functions provide the Power curve equation, that is also known
#' as the Freundlich equation and it is very used in agricultural chemistry,
#' e.g. to model the sorption of xenobiotics in soil. It is also used to model
#' the number of plant species as a function of sampling area
#' (Muller-Dumbois method). These functions provide the equation
#' ('powerCurve.fun()') as well as the self-starters
#' for the \code{\link{nls}} function ( 'NLS.powerCurve()' ) and for the
#' \code{\link[drc]{drm}} function in the 'drc' package ('DRC.powerCurve()')
#'
#'
#' @name SSpowerCurve
#' @aliases powerCurve.fun
#' @aliases NLS.powerCurve
#' @aliases DRC.powerCurve
#'
#' @usage powerCurve.fun(predictor, a, b)
#' NLS.powerCurve(predictor, a, b)
#' DRC.powerCurve(fixed = c(NA, NA), names = c("a", "b"))
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param a model parameter
#' @param b model parameter
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.
#' @param names names. A vector of character strings giving the names of the parameters. The default is reasonable.

#'
#' @details
#' These functions provide the Power curve equation, that is parameterised as:
#' \deqn{ f(x) = a \, x^b }
#' which is totally equivalent to an exponential curve on the logarithm
#' of X:
#' \deqn{ f(x) = a \, \exp \left[ b \, \log(x) \right] }
#' We see that both parameters relate to the ‘slope’ of the curve and b
#' dictates its shape. If 0 < b < 1, the response Y increases as X
#' increases and the curve is convex up. If b < 0 the curve is concave
#' up and Y decreases as X increases. Otherwise, if b > 1, the curve is
#' concave up and Y increases as X increases.
#'
#' @return powerCurve.fun() and NLS.powerCurve() return a numeric value,
#' while DRC.powerCurve() returns a list containing the nonlinear function, the self starter function
#' and the parameter names.
#'
#' @author Andrea Onofri
#'
#' @references Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
#' @references Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
#'
#' @examples
#' dataset <-getAgroData("speciesArea")
#'
#' #nls fit
#' model <- nls(numSpecies ~ NLS.powerCurve(Area, a, b),
#'              data = dataset)
#' summary(model)
#' # drm fit
#' model <- drm(numSpecies ~ Area, fct = DRC.powerCurve(),
#'              data = dataset)
#' summary(model)
#'
#Power Curve ########################################################
# Independently from b, the curve is 0 for x = 0
# The second form adds a displacement on Y axis,
# so that y != 0 when x = 0. Still to be worked
powerCurve.fun <- function(predictor, a, b) {
  a * ( predictor ^ b )
}

powerCurveNO.fun <- function(predictor, a, b, c) {
  a * ( predictor ^ b ) + c
}


powerCurve.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  lmFit <- lm(log(xy[, "y"]) ~ log(xy[, "x"]))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  value <- c(a, b)
  names(value) <- mCall[c("a", "b")]
  value
}

NLS.powerCurve <- selfStart(powerCurve.fun, powerCurve.Init, parameters=c("a", "b"))

# powerCurveNO.Init <- function(mCall, LHS, data, ...) {
#   xy <- sortedXyData(mCall[["predictor"]], LHS, data)
#   pseud
#   pseudoY <- log(xy[, "y"])
#   pseudoX <- log(xy[, "x"])
#   lmFit <- lm(pseudoY ~ pseudoX)
#   coefs <- coef(lmFit)
#   a <- exp(coefs[1])
#   b <- coefs[2]
#   value <- c(a, b)
#   names(value) <- mCall[c("a", "b")]
#   value
# }
#
# NLS.powerCurveNO <- selfStart(powerCurve.fun, powerCurve.Init, parameters=c("a", "b"))


DRC.powerCurve <- function(fixed = c(NA, NA), names = c("a", "b"))
{
  ## Checking arguments
  numParm <- 2
  if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

  ## Fixing parameters (using argument 'fixed')
  notFixed <- is.na(fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]

  ## Defining the non-linear function
  fct <- function(x, parm)
  {


    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm

    a <- parmMat[, 1]; b <- parmMat[, 2]
    a * x ^(b)
  }

  ## Defining self starter function
  ssfct <- function(dataf)
  {
    x <- dataf[, 1]
    y <- dataf[, 2]

    #regression on pseudo y values
    pseudoY <- log( y + 0.00001)
    pseudoX <- log(x)
    coefs <- coef( lm(pseudoY ~ pseudoX) )
    a <- exp(coefs[1])

    b <- coefs[2]

    return(c(a, b)[notFixed])
  }

  ## Defining names
  pnames <- names[notFixed]

  ## Defining derivatives
  deriv1 <- function(x, parms){
      parmMat <- matrix(parmVec, nrow(parms),
                        numParm, byrow = TRUE)
      parmMat[, notFixed] <- parms

      # Approximation by using finite differences
      a <- as.numeric(parmMat[,1])
      b <- as.numeric(parmMat[,2])

      d1.1 <- expoDecay.fun(x, a, b)
      d1.2 <- expoDecay.fun(x, (a + 10e-7), b)
      d1 <- (d1.2 - d1.1)/10e-7

      d2.1 <- expoDecay.fun(x, a, b)
      d2.2 <- expoDecay.fun(x, a, (b + 10e-7) )
      d2 <- (d2.2 - d2.1)/10e-7

      cbind(d1, d2)[notFixed]
    }

    ## Defining the first derivative (in x=dose)
    derivx <- function(x, parm)
    {
      parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
      parmMat[, notFixed] <- parm

      a <- as.numeric(parmMat[,1])
      b <- as.numeric(parmMat[,2])

      d1.1 <- expoGrowth.fun(x, a, b)
      d1.2 <- expoGrowth.fun((x + 10e-7), a, b)
      d1 <- (d1.2 - d1.1)/10e-7
      d1
    }

  ## Defining the ED function

  ## Defining the inverse function

  ## Defining descriptive text
  text <- "Power curve (Freundlich equation)"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text,
                     noParm = sum(is.na(fixed)),
                     deriv1 = deriv1, derivx = derivx)

  class(returnList) <- "drcMean"
  invisible(returnList)
}

"DRC.powerCurveNO" <- function(fixed = c(NA, NA, NA), names = c("a", "b", "c"))
{
  ## Checking arguments
  numParm <- 3
  if (!is.character(names) | !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}

  ## Fixing parameters (using argument 'fixed')
  notFixed <- is.na(fixed)
  parmVec <- rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]

  ## Defining the non-linear function
  fct <- function(x, parm)
  {


    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm

    a <- parmMat[, 1]; b <- parmMat[, 2]; c <- parmMat[, 3]
    a * x ^(b) + c
  }

  ## Defining self starter function
  # ssfct <- function(dataf)
  # {
  #   x <- dataf[, 1]
  #   y <- dataf[, 2]
  #
  #   #regression on pseudo y values
  #   pseudoY <- log( y + 0.00001)
  #   pseudoX <- log(x)
  #   coefs <- coef( lm(pseudoY ~ pseudoX) )
  #   a <- exp(coefs[1])
  #
  #   b <- coefs[2]
  #
  #   return(c(a, b)[notFixed])
  # }

  ## Defining names
  pnames <- names[notFixed]

  ## Defining derivatives

  ## Defining the ED function

  ## Defining the inverse function

  ## Defining descriptive text
  text <- "Power curve not passing for origin"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, names = pnames, text = text, noParm = sum(is.na(fixed)))

  class(returnList) <- "drcMean"
  invisible(returnList)
}


