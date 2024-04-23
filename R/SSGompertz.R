#' Gompertz equations
#'
#' These functions provide the Gompertz equations with 4 (G4.fun), 3 (G3.fun)
#' and 2 (G2.fun) parameters with self-starter for the \code{\link{nls}}
#' function (NLS.G4, NLS.G3 and NLS.G2).
#'
#'
#' @name SSGompertz
#' @aliases G4.fun
#' @aliases G3.fun
#' @aliases G2.fun
#' @aliases NLS.G4
#' @aliases NLS.G3
#' @aliases NLS.G2
#'
#' @usage G4.fun(predictor, b, c, d, e)
#' G3.fun(predictor, b, d, e)
#' G2.fun(predictor, b, e)
#' NLS.G4(predictor, b, c, d, e)
#' NLS.G3(predictor, b, d, e)
#' NLS.G2(predictor, b, e)
#'
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param b model parameter (slope at inflection point)
#' @param c model parameter (lower asymptote)
#' @param d model parameter (higher asymptote)
#' @param e model parameter (abscissa at inflection point)
#'
#' @details
#' The Gompertz equation is parameterised as:
#' \deqn{ f(x) = c + (d - c) \, \exp \left[-exp(-b (x - e))\right] }
#' It is a sygmoidally shaped curve and it is asymmetric about its inflection
#' point. For the 3- and 2-parameter model c is equal to 0, while for the
#' 2-parameter model d is equal to 1.

#' @return  All these functions return a numeric value.
#'
#' @author Andrea Onofri
#'
#' @examples
#' data(beetGrowth)
#' mod3 <- nls(weightInf ~ NLS.G3(DAE, b, c, d), data = beetGrowth)
#' summary(mod3)
#' plot(mod3)

# Gompertz equation for bioassay work
G4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c) * (exp (- exp ( - b * ( x - e))))
}

G4.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- max(y) * 1.05
  c <- min(y) * 0.95

  ## Linear regression on pseudo y values
  pseudoY <- log(-log((y - c)/(d - c)))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b,c,d,e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.G4 <- selfStart(G4.fun, G4.Init, parameters=c("b", "c", "d", "e"))

G3.fun <- function(predictor, b, d, e) {
  x <- predictor
  d * (exp (- exp ( - b * ( x - e))))
}

G3.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  # y <- beetGrowth$weightFree; x <- beetGrowth$DAE
  d <- max(y) * 1.05
  # print(d)
  ## Linear regression on pseudo y values
  pseudoY <- log(-log(y/d))
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- -coefs[2]
  e <- k/b
  value <- c(b,d,e)
  names(value) <- mCall[c("b", "d", "e")]
  value
}

NLS.G3 <- selfStart(G3.fun, G3.Init, parameters=c("b", "d", "e"))

G2.fun <- function(predictor, b, e) {
  x <- predictor
  exp (- exp ( - b * ( x - e)))
}

G2.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]

  ## Linear regression on pseudo y values
  pseudoY <- log(-log(y) )
  coefs <- coef( lm(pseudoY ~ x))
  k <- coefs[1]; b <- - coefs[2]
  e <- k/b
  value <- c(b, e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.G2 <- selfStart(G2.fun, G2.Init, parameters=c("b", "e"))
