#' Bragg's Equation
#'
#' These functions provide the Bragg's equations, that is based on
#' the normal (Gaussian) distribution and it supports a maximum,
#' a minimum and inflection points. These functions provide the
#' equations with 4 (bragg.4.fun) and 3 (bragg.3.fun) parameters
#' with self-starters for the \code{\link{nls}}
#' function (NLS.bragg.4, NLS.bragg.3) and the self-starters for
#' the \code{\link[drc]{drm}} function in the drc package (DRC.bragg.4, DRC.bragg.3)
#'
#'
#' @name SSbragg
#' @aliases bragg.4.fun
#' @aliases bragg.3.fun
#' @aliases NLS.bragg.4
#' @aliases NLS.bragg.3
#' @aliases DRC.bragg.4
#' @aliases DRC.bragg.3
#'
#' @usage bragg.4.fun(X, b, c, d, e)
#' bragg.3.fun(X, b, d, e)
#' NLS.bragg.4(X, b, c, d, e)
#' NLS.bragg.3(X, b, d, e)
#' DRC.bragg.4()
#' DRC.bragg.3()
#'
#' @param X a numeric vector of values at which to evaluate the model
#' @param b model parameter (relates to slope at inflection point)
#' @param c model parameter (lower asymptote)
#' @param d model parameter (maximum value)
#' @param e model parameter (abscissa at maximum value)
#'
#'
#' @details
#' The Bragg's equation is parameterised as:
#' \deqn{ f(x) = c + \left(d - c \right) \, \exp(- b \cdot (X - e)^2) }
#' for the 4-parameter model. For the 3-parameter model, c is equal to 0.
#'  It depicts a bell-shaped curve

#' @return  bragg.4.fun, bragg.3.fun, NLS.bragg.4 and NLS.bragg.3
#' return a numeric value, while DRC.bragg.4 and DRC.bragg.3 return
#' a list containing the nonlinear function and the self starter function
#'
#' @author Andrea Onofri
#'
#' @references Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
#' @references Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
#'
#' @examples
#' library(statforbiology)
#' X <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#' Y1 <- c(0.1, 2, 5.7, 9.3, 19.7, 28.4, 20.3, 6.6, 1.3, 0.1)
#' Y2 <- Y1 + 2
#'
#' # nls fit
#' mod.nls <- nls(Y1 ~ NLS.bragg.3(X, b, d, e) )
#' mod.nls2 <- nls(Y2 ~ NLS.bragg.4(X, b, c, d, e) )
#'
#' # drm fit
#' mod.drc <- drm(Y1 ~ X, fct = DRC.bragg.3() )
#' mod.drc2 <- drm(Y2 ~ X, fct = DRC.bragg.4() )
#' plot(mod.drc, ylim = c(0, 30), log = "")
#' plot(mod.drc2, add = TRUE, col = "red")
#'
#'
# Bragg's equation
bragg.3.fun <- function(X, b, d, e){
  d * exp(- b * (X - e)^2)
}

DRC.bragg.3 <- function(){
  fct <- function(x, parm) {
    bragg.3.fun(x, parm[,1], parm[,2], parm[,3])
  }
  ssfct <- function(data){
    # Get the data
    x <- data[, 1]
    y <- data[, 2]

    d <- max(y)
    e <- x[which.max(y)]

    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( (y + 0.0001) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, d, e)
    return( start )
  }
  names <- c("b", "d", "e")
  text <- "Bragg equation with three parameters"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

bragg.3.init <- function(mCall, LHS, data, ...) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]

    d <- max(y)
    e <- x[which.max(y)]

    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( (y + 0.0001) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, d, e)
    names(start) <- mCall[c("b", "d", "e")]
    start
}

NLS.bragg.3 <- selfStart(bragg.3.fun, bragg.3.init, parameters=c("b", "d", "e"))


bragg.4.fun <- function(X, b, c, d, e){
  c + (d - c) * exp(- b * (X - e)^2)
}

DRC.bragg.4 <- function(){
  fct <- function(x, parm) {
    bragg.4.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])
  }
  ssfct <- function(data){
    # Get the data
    x <- data[, 1]
    y <- data[, 2]

    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]


    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( ((y + 0.0001) - c) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, c, d, e)
    return( start )
  }
  names <- c("b", "c", "d", "e")
  text <- "Bragg equation with four parameters"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

bragg.4.init <- function(mCall, LHS, data, ...) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]

    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]


    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- log( ((y + 0.0001) - c) / d )
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- - coefs[1]
    start <- c(b, c, d, e)
    names(start) <- mCall[c("b", "c", "d", "e")]
    start
}

NLS.bragg.4 <- selfStart(bragg.4.fun, bragg.4.init, parameters=c("b", "c", "d", "e"))
