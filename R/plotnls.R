plotnls <- function(x, type = c("average", "all"),
                     xlim = NULL, gridsize = 100,
                     which = 3, ...){
    fm <- x
    res <- gridsize
    type <- match.arg(type)
    if (!inherits(fm, "nls"))
    stop("use only with \"nls\" objects")

    if (!inherits(eval(fm$data), "data.frame"))
    stop("Can only plot models where variable are passed by setting 'data = data.frame'")

    environment(fm$convInfo)

    if(!is.numeric(which) || any(which < 1) || any(which > 3))
    stop("'which' must be in 1:3")

  if(which == 3){
    if(!inherits(fm, "nlsbc")){
      dframe <- eval(fm$data)
      mm <- fm$m
      cc <- fm$call
      pnms <- names(mm$getPars())
      form <- cc$formula
      rhsnms <- all.vars(form[[3]])
      vnms <- rhsnms[!(rhsnms %in% pnms)]
      if (length(vnms) > 1)
            stop("plotnls not yet implemented for > 1 covariate")
      namList <- list(x = as.name(vnms), y = form[[2]])
      x <- dframe[,as.character(namList$x)]
      y <- dframe[,as.character(namList$y)]
    } else {
      dframe <- eval(fm$data)
      mm <- fm$m
      cc <- fm$oldCall
      pnms <- names(mm$getPars())
      form <- cc$formula
      rhsnms <- all.vars(form[[3]])
      vnms <- rhsnms[!(rhsnms %in% pnms)]
      if (length(vnms) > 2)
            stop("plotnls not yet implemented for > 1 covariate")
        namList <- list(x = as.name(vnms), y = form[[2]])
      x <- dframe[,as.character(namList$x)]
      y <- fm$transy

    }

    if(type == "average"){
                y <- tapply(y, list(factor(x)), mean)
                x <- tapply(x, list(factor(x)), mean)
                }

    if(is.null(xlim)) {
      xmin <- min(x)
      xmax <- max(x)
    } else {
      xmin <- xlim[1]
      xmax <- xlim[2]
    }
    step <- (xmax - xmin)/res
    xseq <- seq(xmin, xmax, step)
    xseqDf <- data.frame(xseq)
    names(xseqDf) <- as.character(namList$x)

    if(!inherits(fm, "nlsbc")){
      newData <- predict(fm, newdata = xseqDf)
      plot(y ~ x, data = dframe, ...)
      points(newData ~ xseq, type = "l", ...)
    } else {
      stop("Not yet implemented for nlsbc objects")
    }


  } else {
    #Get values
    x <- fm
    r <- residuals(x)
    yh <- predict(x) # != fitted() for glm
    w <- weights(x)

    if(!is.null(w)) { # drop obs with zero wt: PR#6640
      wind <- w != 0
      r <- r[wind]
      yh <- yh[wind]
      w <- w[wind]
      labels.id <- labels.id[wind]
    }

    n <- length(r)
    if(which == 2){
      s <- sqrt(deviance(x)/df.residual(x))
      ylab23 <- "Standardized residuals"
      r.w <- if (is.null(w)) r else sqrt(w) * r
      rs <- r.w/s
      }

    l.fit <- "Fitted values"

    ##---------- Do the individual plots : ----------
    if (which == 1) {
      ylim <- range(r, na.rm=TRUE)
      ylim <- extendrange(r = ylim, f = 0.08)
      plot(yh, r, xlab = l.fit, ylab = "Residuals", main = "Residuals vs Expected",
           ylim = ylim, type = "p", ...)
      abline(h = 0, lty = 3, col = "black")
    }
    if (which == 2) { ## Normal
      ylim <- range(rs, na.rm=TRUE)
      ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
      # x.coord <- qnorm(ppoints(y))
      # y.coord <- scale(y, scale = T)
      # plot(y.coord ~ x.coord, main = "Manual QQ-plot",
      # xlab = "Theoretical quantiles", ylab = "Standardised values")
      qq <- qqnorm(rs, main = "Normal Q-Q plot", ylab = ylab23, ylim = ylim, ...)
      qqline(rs, lty = 3, col = "gray50")
    }
  }

}

