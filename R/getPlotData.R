#' Get the data for plotting with ggplot()
#'
#' This method works on a model object and retrieves the data for plotting
#' the observed values (average or all data) and predictions from model fit.
#' It is mainly meant to be used with ggplot()
#'
#' @param obj A fitted model object. Methods are provided for drc, drcte and nls objects
#' @param xlim a vector. The interval for the predictor in which predictions are to be obtained
#' @param gridsize numeric. Number of points in the grid (within xlim) used for predictions.
#' @param type a character string specifying whether all the observed points should be plotted (type = "all")
#'     or whether the response should be averaged over the levels of the predictor
#'     (type = "all"). It is disregarded for drcte objects
#' @param npmle.type a character string. For drcte objects, the NPMLE of the cumulative
#'     density function is only specified at the end of each inspection interval, while
#'      it is not unique within each interval. This argument specifies how the CDF
#'      increases within each interval: possible values are "interpolation" (it is
#'      assumed that the CDF increases progressively), "left" (the CDF increases at
#'      the beginning of each interval), "right" (the CDF increases at the end of each
#'      interval; it is very common in survival analysis) and "midpoint" (the CDF
#'      increases in the middle of each interval; it is very common in survival
#'      analysis). This argument is neglected with nls and drc objects.
#' @param ... Other additional arguments.
#'
#' @return This function returns a list of two elements: 'plotPoints' is a dataframe containing
#'     the observed data, 'plotFits' is a dataframe containing model predictions
#' @author Andrea Onofri
#'
#' @examples
#' fileName <- "https://www.casaonofri.it/_datasets/metamitron.csv"
#' metamitron <- read.csv(fileName)
#' mod <- drm(Conc ~ Time, fct = EXD.2(),
#'            data = metamitron, curveid = Herbicide)
#' retVal <- getPlotData(mod)
#' library(ggplot2)
#' ggplot() +
#'   geom_point(data = retVal$plotPoints, mapping = aes(x = Time, y = Conc)) +
#'   geom_line(data = retVal$plotFits, mapping = aes(x = Time, y = Conc)) +
#'   facet_wrap(~ Herbicide)
#'
getPlotData <- function (obj, ...) UseMethod("getPlotData")

#' @rdname getPlotData
getPlotData.drc <- function(obj, xlim = NULL,
                     gridsize = 100,
                     type = c("average", "all"),
                     ... )
{
  # This function is used to explore 'drc' objects
  # and extract the data to be used for ggplots
  # Date of editing: 29/02/2024
  # x <- mod1

    object <- obj
    type <- match.arg(type)

    ## Determining logarithmic scales
    logX <- FALSE

    ## Constructing the plot data
    dataList <- object[["dataList"]]
    dose <- dataList[["dose"]]
    resp <- dataList[["origResp"]]
    curveid <- dataList[["curveid"]]
    plotid <- dataList[["plotid"]]

    ## Modifying "response" in case of SSD
    if (identical(object[["type"]], "ssd"))
    {
        dose <- unlist(with(dataList, tapply(dose, curveid, function(x){sort(x)}))[unique(dataList[["curveid"]])])
        resp <- unlist(with(dataList, tapply(dose, curveid, function(x){ppoints(x, 0.5)}))[unique(dataList[["curveid"]])])
    }

    ## Slight modifications for event-data (only drc objects)
    if (!is.null(plotid))
    {  # used for event-time data
       assayNoOld <- as.vector(plotid)
    } else {
       assayNoOld <- as.vector(curveid)
    }
    uniAss <- unique(assayNoOld)
    numAss <- length(uniAss)

    doPlot <- TRUE

    plotFct <- (object$"curve")[[1]]
    logDose <- (object$"curve")[[2]]
    naPlot <- ifelse(is.null(object$"curve"$"naPlot"), FALSE, TRUE)

    dlNames <- dataList[["names"]]
    doseName <- dlNames[["dName"]]
    respName <- dlNames[["orName"]]
    levNames <- as.character(unique(curveid))

    ## Determining range of dose values
    if (missing(xlim))
    {
        xLimits <- c(min(dose), max(dose))
    } else {
        xLimits <- xlim  # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
    }

    if ((is.null(logDose)) && (logX))
    {
       dosePts <- exp(seq(log(xLimits[1]), log(xLimits[2]), length = gridsize))
       ## Avoiding that slight imprecision produces dose values outside the dose range
       ## (the model-robust predict method is sensitive to such deviations!)
       dosePts[1] <- xLimits[1]
       dosePts[gridsize] <- xLimits[2]
    } else {
       dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)
    }
#    } else {
    # No handling of multi-dimensional dose values
    # }

    # Handling multiple covariates (to be done)
    doseDim <- 1
    # if(!is.vector(dose)){
    #   doseDim <- length(dose[1,])
    #   doseOld <- dose
    #   addVars <- dose[,-1]
    #   dose <- dose[,1]
    # }

    ## Finding minimum and maximum on response scale
    if (is.null(logDose))
    {
        plotMat <- plotFct(dosePts)
    } else {
        plotMat <- plotFct(logDose^(dosePts))
    }

    maxR <- max(resp)
    maxPM <- apply(plotMat, 2, max, na.rm = TRUE)
    if (max(maxPM) > maxR) {maxPM <- maxPM[which.max(maxPM)]} else {maxPM <- maxR}
    # options(warn=0)

    # if (missing(ylim))
    # {
    #     if (missing(xlim))
    #     {
    #         yLimits <- c(min(resp), maxPM)
    #     } else {
    #         yLimits <- getRange(dose, resp, xLimits)
    #     }
    # } else {
    #     yLimits <- ylim
    # }


    ## Cutting away original x values outside the limits
    eps1 <- 1e-8
    logVec <- !( (dose < xLimits[1] - eps1) | (dose > xLimits[2] + eps1) )
    dose <- dose[logVec]
    resp <- resp[logVec]
    assayNoOld <- assayNoOld[logVec]
    level <- uniAss
    lenlev <- length(level)

    retData <- data.frame(dosePts, as.data.frame(plotMat))
    colnames(retData) <- c(doseName, respName)
    row.names(retData) <- 1:length(retData[,1])

    if(numAss == 1) {
      # Prepare observed points
      plotPoints <- switch(type,
                           "all" = cbind(dose, resp),
                           "average" = cbind(tapply(dose, dose, mean),
                                             tapply(resp, dose, mean))
        )
      plotPoints <- as.data.frame(plotPoints)
      colnames(plotPoints) <- c(doseName, respName)
    } else {
      # With multiple curves melt predictions ...
      colnames(retData) <- c(doseName, levNames)
      retData <- tidyr::pivot_longer(retData, names_to = dlNames$cNames,
                          values_to = "CDF",
                          cols = c((doseDim + 1):length(retData[1,])))
      retData <- as.data.frame(retData)
      colnames(retData)[3] <- respName

      # ... and prepare observed points
      plotPoints <- data.frame(dose, curveid, resp)
      if(type == "average"){
        plotPoints <- aggregate(plotPoints[,3], list(plotPoints[,1], plotPoints[,2]), mean)
        row.names(plotPoints) <- 1:length(plotPoints[,1])
      }
      plotPoints <- as.data.frame(plotPoints)
      colnames(plotPoints) <- c(doseName, dlNames$cNames, respName)
      row.names(plotPoints) <- 1:length(plotPoints[,1])
    }

    returnList <- list(plotPoints = plotPoints, plotFits = retData)
    return(returnList)
}

#' @rdname getPlotData
getPlotData.nls <- function(obj, xlim = NULL,
                            gridsize = 100,
                            type = c("average", "all"),...){
    fm <- obj
    type <- match.arg(type)
    dframe <- data.frame(eval(fm$data))
    mm <- fm$m
    cc <- fm$call
    pnms <- names(mm$getPars())
    form <- cc$formula
    rhsnms <- all.vars(form[[3]])
    vnms <- rhsnms[!(rhsnms %in% pnms)]

    # Define the plot points
    if (length(vnms) > 1){
       stop("This method does not yet work with grouped data")
       namList <- list(x = as.name(vnms[[1]]), y = form[[2]], cov = as.name(vnms[[3]]))
       x <- dframe[,as.character(namList$x)]
       y <- dframe[,as.character(namList$y)]
       cov <- dframe[,as.character(namList$cov)]
       if(type == "average"){
         plotPoints <- aggregate(y, list(factor(x), factor(cov)), mean)
         colnames(plotPoints) <- as.character(namList[c(1,3,2)])
         }
    } else {
      namList <- list(x = as.name(vnms), y = form[[2]])
      x <- dframe[,as.character(namList$x)]
      y <- dframe[,as.character(namList$y)]
      if(type == "average"){
                y <- tapply(y, list(factor(x)), mean)
                x <- tapply(x, list(factor(x)), mean)
                # print(y); print(x)
      }
      plotPoints <- data.frame(x = x, y = y)
      names(plotPoints) <- c(deparse(namList$x), deparse(namList$y))
    }

    # Make predictions
    if(is.null(xlim)) {
      xmin <- min(x)
      xmax <- max(x)
    } else {
      xmin <- xlim[1]
      xmax <- xlim[2]
    }
    step <- (xmax - xmin)/gridsize
    xseq <- seq(xmin, xmax, step)


    if (length(vnms) > 1){
      xseqDf <- expand.grid(xseq, levels(factor(cov)))
      names(xseqDf) <- as.character(namList$x, namList$cov)
    } else {
      xseqDf <- data.frame(xseq)
      names(xseqDf) <- as.character(namList$x)
    }

    newData <- predict(fm, newdata = xseqDf)

    plotFits <- data.frame(x = xseqDf, y = newData)
    names(plotFits)[2] <- deparse(namList$y)

    returnList <- list(plotPoints = plotPoints, plotFits = plotFits)

    return(returnList)
}

#' @rdname getPlotData
getPlotData.drcte <- function(obj, xlim = NULL,
                     gridsize = 100,
                     npmle.type = c("interpolation", "midpoint", "right", "left", "none"),
                     ... ){
  npmle.type <- match.arg(npmle.type)
  drcte::plotData(obj, xlim,
                     gridsize,
                     npmle.type)

}
