#' Beta equation
#'
#' These functions provide the beta equation, a threshold model that was derived
#' from the beta density function and it was adapted to describe phenomena taking place only
#' within a minimum and a maximum threshold value (threshold model), for example
#' to describe the germination rate (GR, i.e. the inverse of germination time)
#' as a function of temperature. These functions provide the beta
#' equation (beta.fun), the self-starters for the \code{\link{nls}}
#' function (NLS.beta) and the self-starters for
#' the \code{\link[drc]{drm}} function in the drc package (DRC.beta)
#'
#'
#' @name SSbeta
#' @aliases beta.fun
#' @aliases NLS.beta
#' @aliases DRC.beta
#'
#' @usage beta.fun(X, b, d, Xb, Xo, Xc)
#' NLS.beta(X, b, d, Xb, Xo, Xc)
#' DRC.beta()
#'
#' @param X a numeric vector of values at which to evaluate the model
#' @param b model parameter
#' @param d model parameter
#' @param Xb model parameter (base threshold level)
#' @param Xo model parameter (optimal threshold level)
#' @param Xc model parameter (ceiling threshold level)

#'
#' @details
#' This equation is parameterised as:
#'
#' \deqn{ f(x) = max\left( d \left\{ \left( \frac{X - Xb}{Xo - Xb} \right) \left( \frac{Xc - X}{Xc - Xo} \right) ^{\frac{Xc - Xo}{Xo - Xb}} \right\}^b , 0 \right)}
#'
#' It depicts a curve that is equal to 0 for X < Xb, grows up to a maximum,
#' that is attained at X = Xo and decreases down to 0, that is attained at
#' X = Xc and mantained for X > Xc.

#' @return  beta.fun, NLS.beta return a numeric value,
#' while DRC.beta returns a list containing the nonlinear function
#' and the self starter function
#'
#' @author Andrea Onofri
#'
#' @references Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
#' @references Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
#'
#' @examples
#' X <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#' Y <- c(0, 0, 0, 7.7, 12.3, 19.7, 22.4, 20.3, 6.6, 0, 0)
#'
#' model <- nls(Y ~ NLS.beta(X, b, d, Xb, Xo, Xc))
#' summary(model)
#' modelb <- drm(Y ~ X, fct = DRC.beta())
#' summary(modelb)
#' plot(modelb, log = "")

beta.fun <- function(X, b, d, Xb, Xo, Xc){
  .expr1 <-  (X - Xb)/(Xo - Xb)
  .expr2 <- (Xc - X)/(Xc - Xo)
  .expr3 <- (Xc - Xo)/(Xo - Xb)
  #ifelse(temp > tb & temp < tc, (.expr2*.expr1^(1/.expr3))^a, 0)
  ifelse(X > Xb & X < Xc, d * (.expr1*.expr2^.expr3)^b, 0)
}

DRC.beta <- function(){

  fct <- function(x, parm) {
      # function code here
    beta.fun(x, parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
  }
  ssfct <- function(data){
     # Self-starting code here
    x <- data[, 1]
    y <- data[, 2]

    d <- max(y)
    Xo <- x[which.max(y)]
    firstidx <- min( which(y !=0) )
    Xb <- ifelse(firstidx == 1,  x[1], (x[firstidx] + x[(firstidx - 1)])/2)
    secidx <- max( which(y !=0) )
    Xc <- ifelse(secidx == length(y),  x[length(x)], (x[secidx] + x[(secidx + 1)])/2)
    c(1, d, Xb, Xo, Xc)

    }
  names <- c("b", "d", "Xb", "Xo", "Xc")
  text <- "Beta function"

  ## Returning the function with self starter and names
  returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
  class(returnList) <- "drcMean"
  invisible(returnList)
}

beta.init <- function(mCall, LHS, data, ...) {
    xy <- sortedXyData(mCall[["X"]], LHS, data)
    x <-  xy[, "x"]; y <- xy[, "y"]

    #Self starting code ##############
    d <- max(y)
    Xo <- x[which.max(y)]
    firstidx <- min( which(y !=0) )
    Xb <- ifelse(firstidx == 1,  x[1], (x[firstidx] + x[(firstidx - 1)])/2)
    secidx <- max( which(y !=0) )
    Xc <- ifelse(secidx == length(y),  x[length(x)], (x[secidx] + x[(secidx + 1)])/2)
    start <- c(1, d, Xb, Xo, Xc)
    names(start) <- mCall[c("b", "d", "Xb", "Xo", "Xc")]
    start
}

NLS.beta <- selfStart(beta.fun, beta.init, parameters=c("b", "d", "Xb", "Xo", "Xc"))

