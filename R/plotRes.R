# This function plots the residuals from a nonlinear regression fit
# with drm(), in the 'drc' package
plotRes <- function(x, which = 1, ...){
    fm <- x
    if (!inherits(fm, "drc"))
    stop("use only with \"drc\" objects")

    # Get values
    res <- residuals(fm)
    # expec <- predict(x) # != fitted() for glm
    expec <- fitted(fm)
    w <- fm$data$weights

    if(!is.null(w)) { # drop obs with zero wt: PR#6640
      wind <- w != 0
      res <- res[wind]
      expec <- expec[wind]
      w <- w[wind]
      # labels.id <- labels.id[wind]
    }

    n <- length(res)
    if(which == 2){
      # scaling is done by way of the residual standard deviation
      s <- sqrt(deviance(fm)/df.residual(fm))
      ylab23 <- "Standardized residuals"
      res.w <- if (is.null(w)) res else sqrt(w) * res
      ress <- res.w/s
      }

    l.fit <- "Fitted values"

    ##---------- Do the individual plots : ----------
    if (which == 1) {
      plot(res ~ expec, xlab = "Expected values", ylab = "Residuals",
           main = "Residuals vs Expected",
           type = "p", ...)
      abline(h = 0, lty = 3, col = "black")
    }
    if (which == 2) {
      ## Normal QQ-plot
      # ylim <- range(ress, na.rm=TRUE)
      # ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
      x.coord <- qnorm(ppoints(res))[order(order(res))] # sorted values
      y.coord <- ress # scale(res, scale = T)
      plot(y.coord ~ x.coord, main = "Normal QQ-plot",
          xlab = "Theoretical quantiles", ylab = "Standardised values", ...)
      # qq <- qqnorm(ress, main = "Normal Q-Q plot", ylab = ylab23)
      qqline(ress, lty = 3, col = "gray50")
    }

}

