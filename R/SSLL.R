#' Log-logistic equation
#'
#' These functions provide the loglogistic equation, that has a
#' symmetric sygmoidal shape over the logarithm of time and it has been used
#' for bioassay work. These functions provide the 4-, 3- and 2-parameter
#' equations (LL4.fun(), LL3.fun() and LL2.fun()) as well as the self-starters
#' for the \code{\link{nls}} function (NLS.LL4(), NLS.LL3() and NLS.LL2() )
#'
#'
#' @name SSLL
#' @aliases LL4.fun
#' @aliases LL3.fun
#' @aliases LL2.fun
#' @aliases NLS.LL4
#' @aliases NLS.LL3
#' @aliases NLS.LL2
#'
#' @usage LL4.fun(predictor, b, c, d, e)
#' LL3.fun(predictor, b, d, e)
#' LL2.fun(predictor, b, e)
#' NLS.LL4(predictor, b, c, d, e)
#' NLS.LL3(predictor, b, d, e)
#' NLS.LL2(predictor, b, e)
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param b model parameter (slope at inflection point)
#' @param c model parameter (lower asymptote)
#' @param d model parameter (higher asymptote)
#' @param e model parameter (abscissa at inflection point)
#'
#' @details
#' These functions provide the log-logistic equation for bioassay work
#' This equation (4-parameters) is parameterised as:
#' \deqn{ f(x) = c + \frac{d - c}{\exp ( 1 + \exp ( - b\,(\log(x) - \log(e))))} }
#' For the 3- and 2-parameters model, c is equal to 0, while for the 2-parameter
#' model d is equal to 1.
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
#' dataset <- getAgroData("brassica")
#'
#' model <- nls(FW ~ NLS.LL4(Dose, b, c, d, e), data = dataset)
#' model <- nls(FW ~ NLS.LL3(Dose, b, d, e), data = dataset)
#' model <- nls(FW/max(FW) ~ NLS.LL2(Dose, b, e), data = dataset)
#' summary(model)
#'
#Log-Logistic Function for bioassay work nlsLL.4
LL4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c+(d-c)/(1+exp(- b*(log(x+0.000001)-log(e))))
}

#NLSLL.4mean <- deriv(~c+(d-c)/(1+exp(b*(log(predictor+0.000001)-log(ED50)))),c("c","d","b","ED50"),function(predictor,c,d,b,ED50){})

LL4.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  c <- min(y) * 0.95
  d <- max(y) * 1.05

  ## Linear regression on pseudo y values
  pseudoY <- log((d-y)/(y-c))
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
  k <- coefs[1]; b <- - coefs[2]
  e <- exp(k/b)
  value <- c(b,c,d,e)
  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.LL4 <- selfStart(LL4.fun, LL4.Init, parameters=c("b", "c", "d", "e"))

# Log-Logistic Function for bioassay work nlsLL.3
# Edited on 07/02/2020
LL3.fun <- function(predictor, b, d, e) {
                      x <- predictor
                      d/(1+exp(-b*(log(x+0.000001)-log(e))))
}

LL3.init <- function(mCall, LHS, data, ...) {
          xy <- sortedXyData(mCall[["predictor"]], LHS, data)
          x <-  xy[, "x"]; y <- xy[, "y"]
          d <- max(y) * 1.05
          ## Linear regression on pseudo y values
          pseudoY <- log((d-y)/(y+0.00001))
          coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
          k <- coefs[1]; b <- - coefs[2]
          e <- exp(k/b)
          value <- c(b,d,e)
          names(value) <- mCall[c("b", "d", "e")]
          value
}

NLS.LL3 <- selfStart(LL3.fun, LL3.init, parameters=c("b", "d", "e"))

# Log-Logistic Function for bioassay work nlsLL.2
# Edited on 07/02/2020
LL2.fun <- function(predictor, b, e) {
  x <- predictor
  1/(1+exp(-b*(log(x+0.000001)-log(e))))
}

LL2.init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]
  d <- 1
  ## Linear regression on pseudo y values
  pseudoY <- log((d-y)/(y+0.00001))
  coefs <- coef( lm(pseudoY ~ log(x+0.000001)))
  k <- coefs[1]; b <- - coefs[2]
  e <- exp(k/b)
  value <- c(b,e)
  names(value) <- mCall[c("b", "e")]
  value
}

NLS.LL2 <- selfStart(LL2.fun, LL2.init, parameters=c("b", "e"))
