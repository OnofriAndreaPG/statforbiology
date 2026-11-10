biplot.AMMIobject <- function(x, biplot = 1, xlim=NULL, ylim=NULL,
                              elabels=NULL, glabels=NULL, quad=FALSE,
                              cexG=0.9, cexE=0.9, xlab=NULL, ylab=NULL,
                              font=1, ...){
    # Edited on 26/5/2023
    obj <- x
    E <- obj$environment_scores
    G <- obj$genotype_scores
    envir.mean <- obj$environment_means
    var.mean <- obj$genotype_means
    overall.mean <- mean(envir.mean)
    if (biplot == 1) {
        plot(1, type = "n",
            xlim = if(is.null(xlim)) range(c(envir.mean, var.mean)) else xlim,
            ylim = if(is.null(ylim)) range(c(E[, 1], G[, 1])) else ylim,
            xlab = if(is.null(xlab)) "VarY" else xlab,
            ylab = if(is.null(ylab)) "PC 1" else ylab,)
        points(envir.mean, E[, 1], "n", col = "red", lwd = 5)
        text(envir.mean, E[, 1],
            labels = if(is.null(elabels)) sub(' +$', '', as.vector(names(envir.mean))) else elabels,
            #adj = c(0.5, 0.5),
            col = "red", font=font,
            cex = cexE)
        points(var.mean, G[, 1], "n", col = "blue", lwd = 5)
        text(var.mean, G[, 1],
            labels = if(is.null(glabels)) sub(' +$', '', as.vector(names(var.mean))) else glabels,
            #adj = c(0.5, 0.5),
            col = "blue", font=font,
            cex = cexG)
        abline(h = 0, lty = 5)
        abline(v = overall.mean, lty=5)
    }
    else {
        plot(1, type = "n",
        xlim = if(is.null(xlim)) range(c(E[, 1], G[, 1])) else xlim,
        ylim = if(is.null(ylim)) range(c(E[, 2], G[, 2])) else ylim,
        xlab = "PC 1", ylab = "PC 2")
        points(E[, 1], E[, 2], "n", col = "red", lwd = 5)
        text(E[, 1], E[, 2],
        	labels = if(is.null(elabels)) sub(' +$', '', as.vector(row.names(E))) else elabels,
        	col = "red", font=2,
        	cex=cexE)
        points(G[, 1], G[, 2], "n", col = "blue", lwd = 5)
        text(G[, 1], G[, 2],
        labels = if(is.null(glabels)) sub(' +$', '', as.vector(row.names(G))) else glabels,
        	#adj = c(0.5, 0.5),
        	col = "blue", font=font,
        	cex=cexG)
        if (quad == TRUE){
        abline(h = 0, lty = 5)
        abline(v = 0, lty=5)
        }
        # if (percDev != FALSE){
        #   expvar <- c(paste(round(obj$mult_Interaction$Perc_of_Total_SS[1]),"%"), paste(round(obj$mult_Interaction$Perc_of_Total_SS[2]),"%"))
        #   arrows(percDev[1], percDev[2], (percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], length=0.15)
        #   text((percDev[1]+max(c(E[, 1], G[, 1]))/7), percDev[2], label=expvar[1], pos=4)
        #   arrows(percDev[1], percDev[2],percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), length=0.15)
        #   text(percDev[1], (percDev[2]+max(c(E[, 1], G[, 1]))/7), label=expvar[2], pos=3)
        # }
    }
}

biplot.GGEobject <- function(x, biplot = 1, xlim=NULL, ylim=NULL,
                              elabels=NULL, glabels=NULL, quad=FALSE,
                              cexG=0.9, cexE=0.9, xlab=NULL, ylab=NULL,
                              font=1, ...){
  obj <- x
  biplot <- 2
  biplot.AMMIobject(obj, biplot, xlim, ylim, elabels, glabels,
                    quad, cexG, cexE, xlab, ylab, font)

}

biplot.polygonView <- function(obj, xlim=NULL, ylim=NULL,
                              elabels=NULL, glabels=NULL, quad=FALSE,
                              cexG=0.9, cexE=0.9, xlab=NULL, ylab=NULL,
                              font=1, ...){
  # It draws a biplot with which-won-where pattern for GGE analysis
  # it should return an error if there is no existing plot
  # Edited on 10/11/2026
  # if(dev.cur() == 1L ) stop("There is no existing plot")
  if(!inherits(obj, "GGEobject")) stop("The polygon view can only be produced with GGE objects")

  # segment <- function(x1, y1, x2, y2, lty=NULL, col = 1){
  # # This function draws a segment joining two points.
  # curve(((x - x1)*(y2 - y1)) / (x2 - x1) + y1,
  #       from=x1, to=x2, lty = if(is.null(lty)) 1 else lty, add=T, col = col)
  # }
  # perp2p <- function(x1, y1, x2, y2, x3, y3, lty=NULL, col = 1){
  # # It draws the perpendicular to a segment from P1(x1,y1) to P2(x2, y2)
  # # passing through the point P3(x3,y3)
  # m <- (y2 - y1)/(x2 - x1)
  # perp.m <- -1/m
  # # fromP <- if(max(x1,x2)>=x3 & min(y1,y2)<=y3) x3 else NULL
  # if(all(c(x1,x2) >= x3)){
  #   fromP <- x3
  #   toP <- NULL
  # } else if(all(c(x1,x2) <= x3)){
  #   fromP <- NULL
  #   toP <- x3
  # } else if (all(c(y1,y2) >= y3)){
  #   fromP <- NULL
  #   toP <- x3
  # } else {
  #   fromP <- NULL
  #   toP <- x3
  # }
  #
  # # toP <- if(min(x1,x2)<=x3 & max(y1,y2)>=y3) x3 else NULL
  # # toP <- if(all(c(x1,x2)<=x3) ) x3 else NULL
  #
  # curve(-1/m*(x - x3) + y3, from = fromP, to=toP,
  #       lty = if(is.null(lty)) 1 else lty, add=T, col = col)
  # }


  # Draws a biplot
  biplot <- 2
  # Find the vertices, by using the 'chull function.
  # Convex hull: the smallest convex polygon that contains all points in a given set
  E <- obj$environment_scores
  G <- obj$genotype_scores
  vertices <- chull(G[,1], G[,2])
  g_hull <- G[vertices, ]
  obj$genotype_scores <- g_hull
  biplot.AMMIobject(obj, biplot, xlim, ylim, elabels, glabels,
                    quad, cexG, cexE, xlab, ylab, font)

  polygon(g_hull, lty = 2)

  # Now I have to find the perpendiculars
  # 1. Completo il poligono
  g_hull_closed <- rbind(g_hull, g_hull[1, ])
  for (i in 1:(nrow(g_hull_closed) - 1)) {
    x1 <- g_hull_closed[i, 1]
    y1 <- g_hull_closed[i, 2]
    x2 <- g_hull_closed[i + 1, 1]
    y2 <- g_hull_closed[i + 1, 2]

    m <- (y2 - y1)/(x2 - x1)
    q <- y1 - m * x1

    perp.m <- -1/m

    x_int <- q / (perp.m - m)
    y_int <- perp.m * x_int
    x_end <- x_int * 30
    y_end <- y_int * 30
    segments(x0 = 0, y0 = 0, x1 = x_end, y1 = y_end, col = "red", lty = 1)
  }
}


