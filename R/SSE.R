#' Modified Gompertz equations
#'
#' These functions provide the modified Gompertz equations with 4 (E4.fun), 3 (E3.fun)
#' and 2 (E2.fun) parameters with self-starter for the \code{\link{nls}}
#' function (NLS.E4, NLS.E3 and NLS.E2) and for the \code{\link[drc]{drm}} function
#' in the 'drc' package (DRC.E4, DRC.E3 and DRC.E2).
#'
#'
#' @name SSE
#' @aliases E4.fun
#' @aliases E3.fun
#' @aliases E2.fun
#' @aliases NLS.E4
#' @aliases NLS.E3
#' @aliases NLS.E2
#' @aliases DRC.E4
#' @aliases DRC.E3
#' @aliases DRC.E2
#'
#' @usage E4.fun(predictor, b, c, d, e)
#' E3.fun(predictor, b, d, e)
#' E2.fun(predictor, b, e)
#' NLS.E4(predictor, b, c, d, e)
#' NLS.E3(predictor, b, d, e)
#' NLS.E2(predictor, b, e)
#' DRC.E4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
#' DRC.E3(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
#' DRC.E2(fixed = c(NA, NA), names = c("b", "e"))
#'
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param b model parameter (slope at inflection point)
#' @param c model parameter (lower asymptote)
#' @param d model parameter (higher asymptote)
#' @param e model parameter (abscissa at inflection point)
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.
#' @param names names. A vector of character strings giving the names of the parameters. The default is reasonable.

#'
#' @details
#' The modified Gompertz equation is parameterised as:
#' \deqn{ f(x) = c + (d - c) \, (1 - \exp \left[-exp( b (x - e))) \right] }
#' It is a sygmoidally shaped curve and it is asymmetric about its inflection
#' point, but the type of asymmetry is different from the Gompertz equation.
#' For the 3- and 2-parameter model c is equal to 0, while for the
#' 2-parameter model d is equal to 1.

#' @return E4.fun, E3.fun, E2.fun, NLS.E4, NLS.E3 and NLS.E2 return a numeric value,
#' while DRC.E4, DRC.E3 and DRC.E2 return a list containing the nonlinear function,
#' the self starter function and the parameter names.
#'
#' @author Andrea Onofri
#'
#' @examples
#' data(beetGrowth)
#' mod3 <- nls(weightInf ~ NLS.E3(DAE, b, c, d), data = beetGrowth)
#' summary(mod3)
#' plot(mod3)
#'
# Modified gompertz equation for bioassay work
E4.fun <- function(predictor, b, c, d, e) {
  x <- predictor
  c + (d - c) * (1 - exp( - exp ( b * ( x - e))))
}

E4.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]

  y <- as.numeric( tapply(y, factor(x), mean) )
  x <- as.numeric( tapply(x, factor(x), mean) )
  mod <- nls(y ~ NLS.L4(x, b, c, d, e))
  value <- as.numeric(coef(mod))

  names(value) <- mCall[c("b", "c", "d", "e")]
  value
}

NLS.E4 <- selfStart(E4.fun, E4.Init, parameters=c("b", "c", "d", "e"))

E3.fun <- function(predictor, b, d, e) {
  x <- predictor
  d * (1 - exp( - exp ( b * ( x - e))))
}

E3.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]

  y <- as.numeric( tapply(y, factor(x), mean) )
  x <- as.numeric( tapply(x, factor(x), mean) )
  mod <- nls(y ~ NLS.L3(x, b, d, e))
  value <- as.numeric(coef(mod))

  names(value) <- mCall[c("b", "d", "e")]
  value
}

NLS.E3 <- selfStart(E3.fun, E3.Init, parameters=c("b", "d", "e"))

E2.fun <- function(predictor, b, e) {
  x <- predictor
  1 - exp( - exp ( b * ( x - e)))
}

E2.Init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]

  ## Linear regression on pseudo y values
  y <- as.numeric( tapply(y, factor(x), mean) )
  x <- as.numeric( tapply(x, factor(x), mean) )
  mod <- nls(y ~ NLS.L2(x, b, e))
  value <- as.numeric(coef(mod))

  names(value) <- mCall[c("b", "e")]
  value
}

NLS.E2 <- selfStart(E2.fun, E2.Init, parameters=c("b", "e"))

"DRC.E4" <-
  function(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
  {
    ## Checking arguments
    numParm <- 4
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

      E4.fun(x, parm[,1], parm[,2], parm[,3], parm[,4])

    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]

      d <- max(y) * 1.05
      c <- min(y) * 0.95

      ## Linear regression on pseudo y values
      pseudoY <- log(-log((d - y)/(d - c)))
      coefs <- coef( lm(pseudoY ~ x))
      k <- coefs[1]; b <- coefs[2]
      e <- -k/b
      value <- c(b,c,d,e)

      return(value[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Modified Gompertz equation (4 parameters)"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
  }

"DRC.E3" <-
  function(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
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

      E3.fun(x, parm[,1], parm[,2], parm[,3])

    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]

      d <- max(y) * 1.05

      ## Linear regression on pseudo y values
      pseudoY <- log(-log((d - y)/d))
      coefs <- coef( lm(pseudoY ~ x))
      k <- coefs[1]; b <- coefs[2]
      e <- -k/b
      value <- c(b,d,e)

      return(value[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Modified Gompertz equation (3 parameters)"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
  }

"DRC.E2" <- function(fixed = c(NA, NA), names = c("b", "e"))
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

      E2.fun(x, parm[,1], parm[,2])

    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]

      ## Linear regression on pseudo y values
      pseudoY <- log(-log(1.01 - y))
      coefs <- coef( lm(pseudoY ~ x))
      k <- coefs[1]; b <- coefs[2]
      e <- -k/b
      value <- c(b, e)

      return(value[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Modified Gompertz equation (2 parameters)"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames, text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
  }
