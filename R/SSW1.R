#' Weibull equation (Type I)
#'
#' These functions provide the Weibull equation (type I), that has an
#' asymmetric sygmoidal shape and it has been used for bioassay work.
#' These functions provide the 4-, 3- and 2-parameter equations
#' (W1.4.fun(), W1.3.fun() and W1.2.fun()) as well as the self-starters
#' for the \code{\link{nls}} function (NLS.W1.4(), NLS.W1.3() and NLS.W1.2()
#'
#'
#' @name SSW1
#' @aliases W1.4.fun
#' @aliases W1.3.fun
#' @aliases W1.2.fun
#' @aliases NLS.W1.4
#' @aliases NLS.W1.3
#' @aliases NLS.W1.2
#'
#' @usage W1.4.fun(predictor, b, c, d, e)
#' W1.3.fun(predictor, b, d, e)
#' W1.2.fun(predictor, b, e)
#' NLS.W1.4(predictor, b, c, d, e)
#' NLS.W1.3(predictor, b, d, e)
#' NLS.W1.2(predictor, b, e)
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param b model parameter (slope at inflection point)
#' @param c model parameter (lower asymptote)
#' @param d model parameter (higher asymptote)
#' @param e model parameter (abscissa at inlection point)
#'
#' @details
#' These functions provide the Weibull (Type I) equation for bioassay work
#' This equation (4-parameters) is parameterised as:
#' \deqn{ f(x) = c + (d - c) \exp ( - \exp ( - b\,(\log(x) - \log(e)))) }
#' For the 3- and 2-parameters model, c is equal to 0, while for the 2-parameter
#'  model d is equal to 1.
#'
#' @return All these functions return a numeric value
#'
#' @author Andrea Onofri
#'
#' @references Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
#' @references Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
#' @references Ritz, C., Jensen, S.M., Gerhard, D., Streibig, J.C., 2019. Dose-response analysis using R, CRC Press. ed. USA.
#'
#' @examples
#' library(statforbiology)
#' dataset <- getAgroData("brassica")
#' model <- nls(FW ~ NLS.W1.4(Dose, b, c, d, e), data = dataset)
#' model.2 <- nls(FW ~ NLS.W1.3(Dose, b, d, e), data = dataset)
#' model.3 <- nls(FW/max(FW) ~ NLS.W1.2(Dose, b, e), data = dataset)
#' summary(model)
#'
#'
# Weibul type 1 Function for bioassay work
# Edited on 07/02/2020
W1.4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c) * exp ( - exp ( - b*(log(x + 0.000001) - log(e))))
}

W1.4.init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]

  d <- max(y) * 1.05
  c <- min(y) * 0.95

  ## Linear regression on pseudo y values
  pseudoY <- log( - log( (y - c) / (d - c) ) )
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)) )

  b <- - coefs[2]
  e <- exp(coefs[1]/b)
  value <- c(b, c, d, e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.W1.4 <- selfStart(W1.4.fun, W1.4.init, parameters=c("b", "c", "d", "e"))


# Weibul type 2 Function for bioassay work nlsW1.3
# Edited on 07/02/2020
W1.3.fun <- function(predictor, b, d, e) {
                      x <- predictor
                      d * exp ( - exp ( - b*(log(x+0.000001)-log(e))))
}

W1.3.init <- function(mCall, LHS, data, ...) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]

          d <- max(y) * 1.05

          ## Linear regression on pseudo y values
          pseudoY <- log( - log( y / d ) )
          coefs <- coef( lm(pseudoY ~ log(x+0.000001)) )

          b <- - coefs[2]
          e <- exp(coefs[1]/b)
          value <- c(b, d, e)
          names(value) <- mCall[c("b", "d", "e")]
          value
}

NLS.W1.3 <- selfStart(W1.3.fun, W1.3.init, parameters=c("b", "d", "e"))

# Weibul type 2 Function for bioassay work nlsW1.3
# Edited on 07/02/2020
W1.2.fun <- function(predictor, b, e) {
  x <- predictor
  exp ( - exp ( - b*(log(x+0.000001)-log(e))))
}

W1.2.init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]

  d <- 1

  ## Linear regression on pseudo y values
  pseudoY <- log( - log( y / d ) )
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)) )

  b <- - coefs[2]
  e <- exp(coefs[1]/b)
  value <- c(b, e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.W1.2 <- selfStart(W1.2.fun, W1.2.init, parameters=c("b", "e"))

