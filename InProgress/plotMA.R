inout.hex <- function(hbin,mincnt)
{
  if(is.null(hbin@cID)) exit("bin object must have a cID slot, \n try re-binning with ID = TRUE")
  tI <- table(hbin@cID)
  which(hbin@cID%in%(names(tI)[tI<mincnt]))
}

plotMAhex <- function (MA, array = 1, xlab = "A", ylab = "M",
                       main = colnames(MA)[array], 
                       xlim = NULL, ylim = NULL, status = NULL,
                       values, pch, col, cex, nbin=40,
                       zero.weights = FALSE,
                       style = "colorscale", legend = 1.2, lcex = 1,
                       minarea = 0.04, maxarea = 0.8, mincnt = 2,
                       maxcnt = NULL, trans = NULL, inv = NULL,
                       colorcut = NULL,
                       border = NULL, density = NULL, pen = NULL,
                       colramp = function(n){ LinGray(n,beg = 90,end = 15) },
                       newpage = TRUE, type = c("p", "l", "n"),
                       xaxt = c("s", "n"), yaxt = c("s", "n"),
                       verbose = getOption("verbose")) 
{
  if(is.null(main))main <- ""
  switch(class(MA),
      marrayRaw = {
        x <- maA(MA[,array])
        y <- maM(MA[,array])
        w <- maW(MA[,array]) 
      },RGList = {
        MA <- MA.RG(MA[, array])
        array <- 1
        x <- MA$A
        y <- MA$M
        w <- MA$w
    }, MAList = {
        x <- as.matrix(MA$A)[, array]
        y <- as.matrix(MA$M)[, array]
        if (is.null(MA$weights)) 
            w <- NULL
        else w <- as.matrix(MA$weights)[, array]
    }, list = {
        if (is.null(MA$A) || is.null(MA$M)) 
            stop("No data to plot")
        x <- as.matrix(MA$A)[, array]
        y <- as.matrix(MA$M)[, array]
        if (is.null(MA$weights)) 
            w <- NULL
        else
            w <- as.matrix(MA$weights)[, array]
    }, MArrayLM = {
        x <- MA$Amean
        y <- as.matrix(MA$coefficients)[, array]
        if (is.null(MA$weights)) 
            w <- NULL
        else 
            w <- as.matrix(MA$weights)[, array]
    }, matrix = {
        narrays <- ncol(MA)
        if (narrays < 2) 
            stop("Need at least two arrays")
        if (narrays > 5) 
            x <- apply(MA, 1, median, na.rm = TRUE)
        else
            x <- rowMeans(MA, na.rm = TRUE)
        y <- MA[, array] - x
        w <- NULL
    }, ExpressionSet = {
        narrays <- ncol(exprs(MA))
        if (narrays < 2) 
            stop("Need at least two arrays")
        if (narrays > 5) 
            x <- apply(exprs(MA), 1, median, na.rm = TRUE)
        else
            x <- rowMeans(exprs(MA), na.rm = TRUE)
        y <- exprs(MA)[, array] - x
        w <- NULL
        if (missing(main)) 
            main <- colnames(exprs(MA))[array]
    }, stop("MA is invalid object"))
    if (!is.null(w) && !zero.weights) {
        i <- is.na(w) | (w <= 0)
        y[i] <- NA
    }
    if (is.null(xlim)) 
        xlim <- range(x, na.rm = TRUE)
    if (is.null(ylim)) 
        ylim <- range(y, na.rm = TRUE)

    hbin <- hexbin(x,y,xbins=nbin,xbnds=xlim,ybnds=ylim, IDs = TRUE)
    hp <- plot(hbin,xlab = xlab, ylab = ylab, main = main,
               type='n',newpage=newpage)
    ## plot the hexagons
    pushHexport(hp$plot.vp)
    if(is.null(maxcnt)) maxcnt <- max(hbin@count)
    if(is.null(colorcut)) colorcut<-seq(0, 1, length = min(17, maxcnt))
    grid.hexagons(hbin, style=style, minarea = minarea, maxarea = maxarea,
             mincnt = mincnt, maxcnt= maxcnt, trans = trans,
             colorcut = colorcut, density = density, border = border,
             pen = pen, colramp = colramp)
    
    if (is.null(status) || all(is.na(status))) {
        if (missing(pch)) 
            pch <- 16
        if (missing(cex)) 
            cex <- 0.3
        if (missing(col)){
          clrs <- colramp(length(colorcut)-1)
          col <- clrs[1]
        }
        pp <- inout.hex(hbin,mincnt)
        grid.points(x[pp], y[pp], pch = pch[[1]],
                    gp=gpar(cex = cex[1], col=col, fill=col))
    }
    else {
        if (missing(values)) {
            if (is.null(attr(status, "values"))) 
                values <- names(sort(table(status), decreasing = TRUE))
            else 
                values <- attr(status, "values")
        }
        sel <- !(status %in% values)
        nonhi <- any(sel)
        if (nonhi) grid.points(x[sel], y[sel], pch = 16, gp=gpar(cex = 0.3))
        nvalues <- length(values)
        if (missing(pch)) {
            if (is.null(attr(status, "pch"))) 
                pch <- rep(16, nvalues)
            else 
                pch <- attr(status, "pch")
        }
        if (missing(cex)) {
            if (is.null(attr(status, "cex"))) {
                cex <- rep(1, nvalues)
                if (!nonhi) 
                  cex[1] <- 0.3
            }
            else cex <- attr(status, "cex")
        }
        if (missing(col)) {
            if (is.null(attr(status, "col")))
                col <- nonhi + 1:nvalues
            else
                col <- attr(status, "col")
        }
        pch <- rep(pch, length = nvalues)
        col <- rep(col, length = nvalues)
        cex <- rep(cex, length = nvalues)
        for (i in 1:nvalues) {
            sel <- status == values[i]
            grid.points(x[sel], y[sel], pch = pch[[i]], gp=gpar(cex = cex[i], 
                col = col[i]))
        }

    }
    popViewport()
    if (legend >0)
      {
        inner <- hexbin:::getPlt(hp$plot.vp, ret.unit="inches", numeric=TRUE)[1]
        inner <- inner/hbin@xbins
        ysize <- hexbin:::getPlt(hp$plot.vp, ret.unit="inches", numeric=TRUE)[2]
        pushViewport(hp$legend.vp)
        grid.hexlegend(legend, ysize=ysize, lcex = lcex, inner = inner,
                       style= style, minarea= minarea, maxarea= maxarea,
                       mincnt= mincnt, maxcnt= maxcnt,
                       trans=trans, inv=inv,
                       colorcut = colorcut,
                       density = density, border = border, pen = pen,
                       colramp = colramp)

            #if (is.list(pch)) 
            #    legend(x = xlim[1], y = ylim[2], legend = values, 
            #      fill = col, col = col, cex = 0.9)
            #else legend(x = xlim[1], y = ylim[2], legend = values, 
            #    pch = pch, , col = col, cex = 0.9)
        popViewport()
      }
    invisible(list(hbin=hbin,plot.vp=hp$plot.vp,legend.vp=hp$legend.vp))
}

hexMA.loess <- function(pMA,span=.4,col='red')
{
  pushHexport(pMA$plot.vp)
  fit<-loess(pMA$hbin@ycm~pMA$hbin@xcm,weights=pMA$hbin@count,span=span)
  grid.lines(seq(0,16,length=200),predict(fit,seq(0,16,length=200)),
             gp=gpar(col=col),default.units='native')
  popViewport()
  invisible(fit)
}

hexMA.abline <- function(pMA, a, b = NULL, h = numeric(0),
                         v = numeric(0), col='black',
                         lty = 1, lwd = 2, ...)
{
  pushHexport(pMA$plot.vp)    
    if (!missing(col)) {
        col.line <- col
    }
    if (!missing(a)) {
        if (inherits(a, "lm")) {
            coeff <- coef(a)
        }
        else if (!is.null(tryCatch(coef(a), error = function(e) NULL)))
            coeff <- coef(a)
        else coeff <- c(a, b)
        if (length(coeff) == 1)
            coeff <- c(0, coeff)
        if (coeff[2] == 0)
            h <- c(h, coeff[1])
        else if (!any(is.null(coeff))) {
            xx <- current.viewport()$xscale
            yy <- current.viewport()$yscale
            x <- numeric(0)
            y <- numeric(0)
            ll <- function(i, j, k, l) (yy[j] - coeff[1] - coeff[2] *
                xx[i]) * (yy[l] - coeff[1] - coeff[2] * xx[k])
            if (ll(1, 1, 2, 1) <= 0) {
                y <- c(y, yy[1])
                x <- c(x, (yy[1] - coeff[1])/coeff[2])
            }
            if (ll(2, 1, 2, 2) <= 0) {
                x <- c(x, xx[2])
                y <- c(y, coeff[1] + coeff[2] * xx[2])
            }
            if (ll(2, 2, 1, 2) <= 0) {
                y <- c(y, yy[2])
                x <- c(x, (yy[2] - coeff[1])/coeff[2])
            }
            if (ll(1, 2, 1, 1) <= 0) {
                x <- c(x, xx[1])
                y <- c(y, coeff[1] + coeff[2] * xx[1])
            }
            if (length(x) > 0)
                grid.lines(x = x, y = y, default.units = "native",
                  gp = gpar(col = col.line, lty = lty, lwd = lwd))
        }
    }
    h <- as.numeric(h)
    v <- as.numeric(v)
    for (i in seq(along = h))
      grid.lines(y = rep(h[i], 2), default.units = "native",
                 gp = gpar(col = col.line, lty = lty, lwd = lwd))
    for (i in seq(along = v))
      grid.lines(x = rep(v[i], 2), default.units = "native",
                 gp = gpar(col = col.line, lty = lty, lwd = lwd))
    popViewport()
}
