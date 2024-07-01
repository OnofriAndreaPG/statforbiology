#' Lorentz equation
#'
#' These functions provide the Lorentz equation with 3 and 4 parameters ('lorentz.3.fun()'
#' and 'lorentz.4.fun()' ), as well as the self-starters for the \code{\link{nls}}
#' function ( 'NLS.lorentz.3()' and 'NLS.lorentz.4()') and for the
#' \code{\link[drc]{drm}} function in the 'drc' package ('DRC.lorentz.3()' and 'DRC.lorentz.4()')
#'
#'
#' @name SSlorentz
#' @aliases lorentz.3.fun
#' @aliases lorentz.4.fun
#' @aliases NLS.lorentz.3
#' @aliases NLS.lorentz.4
#' @aliases DRC.lorentz.3
#' @aliases DRC.lorentz.4

#'
#' @usage lorentz.3.fun(X, b, d, e)
#' lorentz.4.fun(X, b, c, d, e)
#' NLS.lorentz.3(X, b, d, e)
#' NLS.lorentz.4(X, b, c, d, e)
#' DRC.lorentz.3()
#' DRC.lorentz.4()
#'
#' @param X a numeric vector of values at which to evaluate the model
#' @param b model parameter
#' @param c model parameter
#' @param d model parameter
#' @param e model parameter
#'
#' @details
#' These functions provide the Lorentz equation, that is a bell-shaped
#' curve similar to a gaussian density function. It is parameterised as:
#'
#' \deqn{ f(x) = c + \frac{d - c} {( 1 + b \, (X - e)^2) } }
#'
#' The parameter 'e' represents the abscissa of the maximum value,
#' while c is the minimum (asymptotic) response value and d is the maximum
#' response value. The parameter 'b' relates to the slope at inflection
#' point. For the 3-parameters curve, c is equal to 0.
#'
#' @return lorentz.3.fun(), lorentz.4.fun(), NLS.lorentz.3() and NLS.lorentz.4() return a numeric value,
#' while DRC.lorentz.3() and DRC.lorentz.4()  returns a list containing
#' the nonlinear function, the self starter function and the parameter names.
#'
#' @author Andrea Onofri
#'
#' @references Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
#' @references Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
#'
#' @examples
#'
#' X <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#' Y1 <- c(0.1, 2, 5.7, 9.3, 19.7, 28.4, 20.3, 6.6, 1.3, 0.1)
#' Y2 <- Y1 + 2
#'
#' # nls fit
#' mod.nls <- nls(Y1 ~ NLS.lorentz.3(X, b, d, e) )
#' mod.nls2 <- nls(Y2 ~ NLS.lorentz.4(X, b, c, d, e) )
#'
#' # drm fit
#' mod.drc <- drm(Y1 ~ X, fct = DRC.lorentz.3() )
#' mod.drc2 <- drm(Y2 ~ X, fct = DRC.lorentz.4() )
#' plot(mod.drc, ylim = c(0, 30), log = "")
#' plot(mod.drc2, add = TRUE, col = "red")
#'
#'
# Lorentz's equation
lorentz.3.fun <- function(X, b, d, e){
  d / ( 1 + b * (X - e)^2)
}

DRC.lorentz.3 <- function(){
  fct <- function(x, parm) {
    lorentz.3.fun(x, parm[,1], parm[,2], parm[,3])
  }
  ssfct <- function(data){
    # Get the data
    x <- data[, 1]
    y <- data[, 2]

    d <- max(y)
    e <- x[which.max(y)]

    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- ( d - y )/ y
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, d, e)
    return( start )
  }
  names <- c("b", "d", "e")
  text <- "Lorentz equation with three parameters"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

lorentz.3.init <- function(mCall, LHS, data, ...) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]

    d <- max(y)
    e <- x[which.max(y)]

    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- ( d - y )/ y
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, d, e)
    names(start) <- mCall[c("b", "d", "e")]
    start
}

NLS.lorentz.3 <- selfStart(lorentz.3.fun, lorentz.3.init, parameters=c("b", "d", "e"))


lorentz.4.fun <- function(X, b, c, d, e){
  c + (d - c) / ( 1 + b * (X - e)^2)
}

DRC.lorentz.4 <- function(){
  fct <- function(x, parm) {
    lorentz.4.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])
  }
  ssfct <- function(data){
    # Get the data
    x <- data[, 1]
    y <- data[, 2]

    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]


    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- (d - y)/(y - c)
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, c, d, e)
    return( start )
  }
  names <- c("b", "c", "d", "e")
  text <- "Lorentz equation with four parameters"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

lorentz.4.init <- function(mCall, LHS, data, ...) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]

    d <- max(y)
    c <- min(y) * 0.95
    e <- x[which.max(y)]


    ## Linear regression on pseudo-y and pseudo-x
    pseudoY <- (d - y)/(y - c)
    pseudoX <- (x - e)^2
    coefs <- coef( lm(pseudoY ~ pseudoX - 1) )
    b <- coefs[1]
    start <- c(b, c, d, e)
    names(start) <- mCall[c("b", "c", "d", "e")]
    start
}

NLS.lorentz.4 <- selfStart(lorentz.4.fun, lorentz.4.init, parameters=c("b", "c", "d", "e"))
