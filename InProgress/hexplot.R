## --- Grid version !
plot.hexbin <-
    function(x, style = c("colorscale", "centroids", "lattice",
		"nested.lattice", "nested.centroids"),
	     legend = 1, lcex = 1,
	     minarea = 0.04, maxarea = 0.8, mincnt = 1, maxcnt = max(x$cnt),
	     trans = NULL, inv = NULL,
	     colorcut = seq(0, 1, length = min(17, maxcnt)),
	     border = NULL, density = NULL, pen = NULL,
	     colramp = function(n){ LinGray(n,beg = 90,end = 15) },
	     xlab = "", ylab = "", verbose = getOption("verbose"), ...)
{

    ##NL: Added colramp argument to use fancy ramps
    if(!inherits(x,"hexbin")) stop("first argument must be a hexbin object")
    style <- match.arg(style) # check; user can abbreviate
    if(minarea < 0)
	stop("Minimum area must be non-negative")
    if(maxarea > 1)
	warning("Maximum area should be <= 1 this leads to overlapping hexagons")
    if(minarea > maxarea)
	stop("Minarea must be <= maxarea")
    if(length(colorcut) > 1) {
	if(colorcut[1] != 0)
	    stop("Colorcut lower boundary must be 0")
	if(colorcut[length(colorcut)] != 1)
	    stop("Colorcut upper boundary must be 1")
    }
    else {
	colorcut <- 1
    }

    shape <- x$shape
    ## here is where the plotting starts.
    ## plot.new()
    grid.newpage()


    vp.legend <- viewport(x = 0.5, y = 1, w = 1, h = unit(1, "lines"),
                   just = c("centre", "top"))
    push.viewport(vp.legend)




    din <- oldpar$din
    pin <- oldpar$pin
    mai <- oldpar$mai
    fig <- oldpar$fig
    xsize <- pin[1]
    ysize <- pin[2]
    if(is.logical(legend)) {
	if(legend)
	    stop("Give the legend width")
	else legend <- 0
    }
    xsize <- xsize - legend
    if(xsize < 1)
	stop("plot width too small")
    if(shape * xsize <= ysize) {
	center <- (ysize - shape * xsize)/2
	mai[1] <- mai[1] + center
	mai[3] <- mai[3] + center
	mai[4] <- mai[4] + legend
	ysize <- shape * xsize
    }
    else {
	center <- (xsize - ysize/shape)/2
	mai[2] <- mai[2] + center
	mai[4] <- mai[4] + center + legend
	xsize <- ysize/shape
    }
    vp.main <- viewport(x = 0.5, y = 1, w = 1, h = unit(1, "lines"),
                   just = c("centre", "top"))
    push.viewport(vp.main)
    if(legend > 0) {
	if(!is.null(trans) && is.null(inv))
	    stop("Must supply the inverse transformation")
	if(verbose)
	    cat("plot.hexbin( legend > 0):  ... hex.legend()\n")

	#screen(1,new = FALSE)## (MM: not in bin2d)
	mai[2] <- mai[2] + xsize
	mai[4] <- mai[4] - legend
	par(fig = fig, mai = mai, new=TRUE)
	plot.new()
	plot.window(c(0, legend), c(0, ysize))
	     ##type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
	inner <- xsize/x$xbins
	outer <- (inner * sqrt(3))/2
	##switch(style,
	##	 lattice = ,
	##	 centroids = {
	##	     if(length(colorcut) * outer > ysize - 1) {
	##		 warning("Colorcut is being shortened")
	##		 colorcut <- seq(0, 1,
	##				 max(1, floor((ysize - 1)/outer)))
	##	     }
	##	 }
	##	 )
	hex.legend(legend, ysize=ysize, lcex = lcex, inner = inner,
		   style= style, minarea= minarea, maxarea= maxarea,
		   mincnt= mincnt, maxcnt= maxcnt, trans=trans, inv=inv,
		   colorcut = colorcut,
		   density = density, border = border, pen = pen,
		   colramp = colramp)

	mai[2] <- mai[2] - xsize
	mai[4] <- mai[4] + legend
	## MM: bin2d does   par(new = TRUE)

	## }
	## screen(1,new = FALSE)## (MM: not in bin2d)
	## if(legend>0) {

	if(verbose) cat("plot.hexbin(): calling par(new=TRUE,..); plot.new()")

      par(new=TRUE, mai = mai, fig = fig)
      plot.new()
    }

    plot(x$xbnds, x$ybnds, type = "n", xlab = xlab, ylab = ylab, ...)
    hexagons(x, style = style,
	     minarea = minarea, maxarea = maxarea,
	     mincnt = mincnt, maxcnt = maxcnt,
	     trans = trans, colorcut = colorcut, density = density,
	     border = border, pen = pen, colramp = colramp, verbose = verbose)

    invisible(par(no.readonly=TRUE))
} ## plot.hexbin()
