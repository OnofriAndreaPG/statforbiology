#' Logistic equations
#'
#' These functions provide the logistic equations with 4 (L4.fun), 3 (L3.fun)
#' and 2 (L2.fun) parameters with self-starter for the \code{\link{nls}}
#' function (NLS.L4, NLS.L3 and NLS.L2) and the self-starter for logistic
#' function with two parameters for the \code{\link[drc]{drm}} function in the
#' drc package (DRC.L2).
#'
#'
#' @name SSL
#' @aliases L4.fun
#' @aliases L3.fun
#' @aliases L2.fun
#' @aliases NLS.L4
#' @aliases NLS.L3
#' @aliases NLS.L2
#' @aliases DRC.L2
#'
#' @usage L4.fun(predictor, b, c, d, e)
#' L3.fun(predictor, b, d, e)
#' L2.fun(predictor, b, e)
#' NLS.L4(predictor, b, c, d, e)
#' NLS.L3(predictor, b, d, e)
#' NLS.L2(predictor, b, e)
#' DRC.L2(upper = 1, fixed = c(NA, NA), names = c("b", "e"))
#'
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param b model parameter (slope at inflection point)
#' @param c model parameter (lower asymptote)
#' @param d model parameter (higher asymptote)
#' @param e model parameter (abscissa at inflection point)
#' @param upper numeric. For L.2, a upper asymptote different from 1 can be specified.
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.
#' @param names names. A vector of character strings giving the names of the parameters. The default is reasonable.
#'
#'
#' @details
#' The logistic equation is parameterised as:
#' \deqn{ f(x) = c + \frac{d - c}{1+exp\left[-b(x - e)\right]} }
#' for the 3- and 2-parameter model c is equal to 0, while for the 2-parameter model
#' d is equal to 1.

#' @return  L4.fun, L3.fun, L2.fun, NLS.L4, NLS.L3 and NLS.L2 return a numeric value,
#' while DRC.L2 returns a list containing the nonlinear function, the self starter function
#' and the parameter names.
#'
#' @author Andrea Onofri
#'
#' @examples
#' data(beetGrowth)
#' mod3 <- nls(weightInf ~ NLS.L3(DAE, b, c, d), data = beetGrowth)
#' mod3b <- drm(weightInf ~ DAE, fct=DRC.L2(upper = 25), data = beetGrowth)


# Logistic Function for bioassay work nlsL.4
# Similar to L.4() in drc
L4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c)/(1 + exp( - b* (x - e)))
}

L4.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- max(y) * 1.05
  c <- min(y)

  ## Linear regression on pseudo y values
  pseudoY <- log((d - y)/(y+0.00001 - c))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b,c,d,e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.L4 <- selfStart(L4.fun, L4.Init, parameters=c("b", "c", "d", "e"))

# Logistic Function for bioassay work nlsL.3
# Similar to L.3() in drc
L3.fun <- function(predictor, b, d, e) {
  x <- predictor
  d/(1 + exp( - b* (x - e)))
}

L3.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  # y <- beetGrowth$weightFree; x <- beetGrowth$DAE
  d <- max(y) * 1.05
  # print(d)
  ## Linear regression on pseudo y values
  pseudoY <- log((d-y)/(y+0.00001))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- -coefs[2]
  e <- k/b
  value <- c(b,d,e)
  names(value) <- mCall[c("b", "d", "e")]
  value
}

NLS.L3 <- selfStart(L3.fun, L3.Init, parameters=c("b", "d", "e"))

# Logistic Function for bioassay work nlsL.2
# Similar to L.2() in drc
L2.fun <- function(predictor, b, e) {
  x <- predictor
  1/(1 + exp( - b* (x - e)))
}

L2.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- 1
  ## Linear regression on pseudo y values
  pseudoY <- log((d - y)/(y+0.00001))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b, e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.L2 <- selfStart(L2.fun, L2.Init, parameters=c("b", "e"))

# Logistic curve with two-parameters (DRC) ##########################
"DRC.L2" <-
  function (upper = 1, fixed = c(NA, NA), names = c("b", "e"))
  {
    numParm <- 2
    if (!is.character(names) | !(length(names) == numParm)) {
      stop("Not correct 'names' argument")
    }
    if (!(length(fixed) == numParm)) {
      stop("Not correct length of 'fixed' argument")
    }
    return(logistic(fixed = c(fixed[1], 0, upper, fixed[2], 1),
                    names = c(names[1], "c", "d", names[2], "f"),
                    fctName = as.character(match.call()[[1]]),
                    fctText = "Logistic (ED50 as parameter)"))
  }

