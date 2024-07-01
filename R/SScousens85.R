#' Rectangular hyperbola for yield/weed density relationship
#'
#' These functions provide the rectangula hyperbola that was devided by Cousens (1985)
#' for modelling the relationship between crop yield and weed density. The function was
#' derived from the yield-loss function, and contains parameters that are revelant for
#' competition studies. These functions provide the
#' equation (cousens85.fun), the self-starters for the \code{\link{nls}}
#' function (NLS.cousens85) and the self-starters for
#' the \code{\link[drc]{drm}} function in the drc package (DRC.cousens85)
#'
#'
#' @name SScousens85
#' @aliases cousens85.fun
#' @aliases NLS.cousens85
#' @aliases DRC.cousens85
#'
#' @usage cousens85.fun(predictor, Ywf, i, A)
#' NLS.cousens85(predictor, Ywf, i, A)
#' DRC.cousens85(fixed = c(NA, NA, NA), names = c("Ywf", "i", "A"))
#'
#' @param predictor a numeric vector of values at which to evaluate the model
#' @param Ywf model parameter (Weed-free yield)
#' @param i model parameter (initial slope)
#' @param A model parameter (maximum percentage yield loss)
#' @param fixed numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.
#' @param names a vector of character strings giving the names of the parameters. The default is usually reasonable.

#'
#' @details
#' This equation is parameterised as:
#' \deqn{ f(x) = Ywf \, \frac{(1 - (i predictor)} {(100 \, (1 + i \, predictor/A)))} }
#' It depicts a decreasing curve with no inflection point. The curve is equal
#' to 'Ywf' when x = 0 and the lower asymptote is at 'A' multiplied by 'Ywf/100'

#' @return  cousens85.fun, NLS.cousens85 return a numeric value,
#' while DRC.cousens85 return a list containing the nonlinear function
#' and the self starter function
#'
#' @author Andrea Onofri
#'
#' @references Ratkowsky, DA (1990) Handbook of nonlinear regression models. New York (USA): Marcel Dekker Inc.
#' @references Onofri, A. (2020). A collection of self-starters for nonlinear regression in R. See: \url{https://www.statforbiology.com/2020/stat_nls_usefulfunctions/}
#' @references Cousens, R., 1985. A simple model relating yield loss to weed density. Annals of Applied Biology 107, 239–252. https://doi.org/10.1111/j.1744-7348.1985.tb01567.x
#'
#' @examples
#' library(statforbiology)
#' dataset <- getAgroData("Sinapis")
#'
#' # nls fit
#' mod.nls <- nls(yield ~ NLS.cousens85(density, Ywf, i, A),
#'                data = dataset )
#' summary(mod.nls)
#' mod.nls2 <- drm(yield ~ density, fct = DRC.cousens85(), data = dataset )
#' summary(mod.nls2)
#' plot(mod.nls2)
#'
# Yield - Weed Density function ###################
cousens85.fun <- function(predictor, Ywf, i, A) {
  Ywf * (1 - (i * predictor) / (100 * (1 + i * predictor/A)))
}
cousens85.init <- function(mCall, LHS, data, ...) {
  xy <- sortedXyData(mCall[["predictor"]], LHS, data)
  x <- xy[, "x"]; y <- xy[, "y"]

  Ywf <- max(y)+10e-06
  YL <- (1 - y/Ywf)*100
  #regression on pseudo y values
  pseudoY <- 1 /  (YL + 0.000001)
  pseudoX <- 1 / (x + 0.00001)
  coefs <- coef( lm(pseudoY ~ pseudoX) )
  A <- 1 / coefs[1]
  i <- 1 / coefs[2]
  value <- c(Ywf, i, A)
  names(value) <- mCall[c("Ywf", "i", "A")]
  value
}

NLS.cousens85 <- selfStart(cousens85.fun, cousens85.init,
                           parameters=c("Ywf", "i", "A"))

"DRC.cousens85" <-
  function(fixed = c(NA, NA, NA), names = c("Ywf", "i", "A"))
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

      Ywf <- parmMat[, 1]; i <- parmMat[, 2]; A <- parmMat[, 3]
      Ywf*(1 - (i*x) / (100 * (1 + i * x/A)))
    }

    ## Defining self starter function
    ssfct <- function(dataf)
    {
      x <- dataf[, 1]
      y <- dataf[, 2]

      Ywf <- max(y)+10e-06
      YL <- (1 - y/Ywf)*100
      #regression on pseudo y values
      pseudoY <- 1 /  (YL + 0.000001)
      pseudoX <- 1 / (x + 0.00001)
      coefs <- coef( lm(pseudoY ~ pseudoX) )
      A <- 1 / coefs[1]; i <- 1 / coefs[2]

      return(c(Ywf, i, A)[notFixed])
    }

    ## Defining names
    pnames <- names[notFixed]

    ## Defining derivatives

    ## Defining the ED function

    ## Defining the inverse function

    ## Defining descriptive text
    text <- "Yield-Weed Density function (Cousens, 1985)"

    ## Returning the function with self starter and names
    returnList <- list(fct = fct, ssfct = ssfct, names = pnames,
                       text = text, noParm = sum(is.na(fixed)))

    class(returnList) <- "drcMean"
    invisible(returnList)
  }
