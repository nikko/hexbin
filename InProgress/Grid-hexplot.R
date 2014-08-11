## This is an attemp to plot a hexbin object so that if the plot is 
## resized the shapes of the hexagons will be preserved. It may be 
## that this is not possible, so that this will have to be redone to
## accept xy coordinates and bin on the fly so that the shape parameter
## of the bins can vary according to the aspect ratio of the plot.
## This is not desireable since with a huge data set one would ideally
## Like to only bin once.

gplot.hexbin <-
    function(x, style = c("colorscale", "centroids", "lattice",
		"nested.lattice", "nested.centroids"),
	     legend = 1.2, lcex = 1,
	     minarea = 0.04, maxarea = 0.8, mincnt = 1, maxcnt = max(x@count),
	     trans = NULL, inv = NULL, 
	     colorcut = seq(0, 1, length = min(17, maxcnt)),
	     border = NULL, density = NULL, pen = NULL,
	     colramp = function(n){ LinGray(n,beg = 90,end = 15) },
	     xlab = "", ylab = "", main="", newpage=TRUE, 
             verbose = getOption("verbose"), ...)
{

    
    if(!is(x,"hexbin")) stop("first argument must be a hexbin object")

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
        if(colorcut > 1) colorcut <- seq(0, 1,
                                         length = min(c(17,colorcut,maxcnt)))
	else colorcut <- 1
    }

    shape <- x@shape
    ## here is where the plotting starts.
    if(is.logical(legend)) {
	if(legend)
	    stop("Give the legend width")
	else legend <- 0
    }
    if (newpage) grid.newpage()
    
    #din <- c(convertWidth(unit(1,"npc"),"inches"),
             #convertHeight(unit(1,"npc"),"inches"))
    vpin <- c(convertWidth(unit(1,"npc"),"inches"),
           convertHeight(unit(1,"npc"),"inches"))

    fig <- c(vpin[1] - legend,vpin[2])

    ## Now set up the margins and the plotting region of the plot
    ## viewport so that the aspect ratio is correct
    mai <- c(convertY(unit(5.1,"lines"),"inches"),
             convertX(unit(4.1,"lines"),"inches"),
             convertY(unit(4.1,"lines"),"inches"),
             convertX(unit(2.1,"lines"),"inches"))
    
    pin <- c(fig[1]-(mai[2]+mai[4]), fig[2]-(mai[1]+mai[3]))
    xsize <- pin[1]
    ysize <- pin[2]
##
    if(shape * xsize <= ysize) {
	center <- (ysize - shape * xsize)/2
	mai[1] <- mai[1] + center
	mai[3] <- mai[3] + center
	ysize <- shape * xsize
    }
    else {
	center <- (xsize - ysize/shape)/2
	mai[2] <- mai[2] + center
	mai[4] <- mai[4] + center
	xsize <- ysize/shape
    }
    fig <- c(pin[1]+mai[2]+ mai[4],fig[2])
    mar <- c(convertY(unit(mai[1],"inches"),"lines"),
             convertX(unit(mai[2],"inches"),"lines"),
             convertY(unit(mai[3],"inches"),"lines"),
             convertX(unit(mai[4],"inches"),"lines"))

    fig.vp <-viewport(x=0,y=0,height=1,
                        width=convertUnit(unit(fig[1],"inches"),"npc"),
                        just=c("left","bottom")) 
    pushViewport(fig.vp)
    plot.vp <- plotViewport(margins=mar,
                            xscale=x@xbnds + c(-0.05, 0.05) * diff(x@xbnds),
                            yscale=x@ybnds + c(-0.05, 0.05) * diff(x@ybnds)) 
    pushViewport(plot.vp)
    grid.rect()
    grid.xaxis()
    grid.yaxis()
    grid.hexagons(x, style = style,
	     minarea = minarea, maxarea = maxarea,
	     mincnt = mincnt, maxcnt = maxcnt,
	     trans = trans, colorcut = colorcut, density = density,
	     border = border, pen = pen, colramp = colramp, verbose = verbose)
    if(nchar(xlab)>0)
      grid.text(xlab, y = unit(-2, "lines"), gp = gpar(fontsize = 16))
    if(nchar(ylab)>0)
      grid.text(ylab, x = unit(-2, "lines"), gp = gpar(fontsize = 16), rot = 90)
    if(nchar(main)>0)
      grid.text(main, y = unit(1, "npc") + unit(1.5, "lines"),
                gp = gpar(fontsize = 18))

    popViewport()
    popViewport()
    if(legend > 0) {
	if(!is.null(trans) && is.null(inv))
	    stop("Must supply the inverse transformation")
	if(verbose)
	    cat("plot.hexbin( legend > 0):  ... hex.legend()\n")
	inner <- xsize/x@xbins
	##outer <- (inner * sqrt(3))/2
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
        vp.legend <- viewport(x=unit(1,"npc") - 
                                convertX(unit(legend,"inches"),"npc"),
                              y=convertY(unit(mai[1],"inches"),"npc"),
                              height=unit(1,"npc") -
                                convertY(unit(mai[3]+mai[1],"inches"),"npc"),
                              width=convertUnit(unit(legend,"inches"),"npc"),
                              default.units = "native",
                              just=c("left","bottom"),xscale = c(0, legend),
                              yscale=c(0, ysize))
        pushViewport(vp.legend)
	grid.hexlegend(legend, ysize=ysize, lcex = lcex, inner = inner,
		   style= style, minarea= minarea, maxarea= maxarea,
		   mincnt= mincnt, maxcnt= maxcnt, trans=trans, inv=inv,
		   colorcut = colorcut,
		   density = density, border = border, pen = pen,
		   colramp = colramp)

    }
  invisible()
} ## plot.hexbin()
