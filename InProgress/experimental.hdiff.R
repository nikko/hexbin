hull2poly <- function(ind,xy)
{
  m <- matrix(c(xy$x[ind], xy$y[ind]),nc=2) # FIXME: cbind()!
  as(m, "gpc.poly")
}

inbb <- function(ply,pts)
{
  bb <- get.bbox(ply)
  (pts$x > bb$x[1]) & (pts$x < bb$x[2]) & (pts$y > bb$y[1]) & (pts$y < bb$y[2])
}

which.in.poly <- function(pts,poly)
{
  inset <- inbb(poly,pts)
  ppts <- get.pts(poly)
  nvert <- as.integer(length(ppts[[1]]$x)+1)
  npt <- as.integer(length(pts$x))
  #answ<- .C("ppoly", nv=nvert, np=npt,
  #          px=as.double(c(ppts[[1]]$x,ppts[[1]]$x[1])),
  #          py=as.double(c(ppts[[1]]$y,ppts[[1]]$y[1])),
  #          x=as.double(pts$x),y=as.double(pts$y),
  #          inp=as.integer(as.numeric(inset)),package="hexbin")
  #as.logical(answ$inp)
  inset[inset]<- point.in.polygon(pts$x[inset], pts$y[inset],
                                  ppts[[1]]$x, ppts[[1]]$y)
  inset
}

hdiffplot <-
    function(bin1, bin2,
             xbnds1=NULL, ybnds1=NULL,
             xbnds2=NULL, ybnds2=NULL,
             density = rep.int(-1,6), border = rep.int(FALSE, 6), pen = 2:7,
             size = 0.1, lwd = 2, eps = 1e-6, unzoom = 1.04,
             xlab = "", ylab = "", main="", ...)
{
#### FIXME -------------# This should be in ../DESCRIPTION
    require(gpclib)     # 'Depends: ...'
#### FIXME -------------# ***once*** the current file is in the package!!

    ## Arguments:
    ##  bin1    : hexagon bin object
    ##  bin2    : hexagon bin object
    ##            Note: bin objects must have overlapping plotting bounds
    ##                  and the same shape parameter.
    ##  density : for hexagon graph paper
    ##  border  : plot the border of the hexagon, use TRUE for
    ##            hexagon graph paper
    #require(gstat)
    if(bin1@shape!=bin2@shape)
      stop("bin objects must have same shape parameter")
    ##_______________ Collect computing constants______________
    tmp1 <- hcell2xy(bin1)
    tmp2 <- hcell2xy(bin2)
    shape <- bin1$shape
    if(all(bin1@xbnds==bin2@xbnds) & all(bin1@ybnds==bin2@ybnds))
      equal.bounds <- TRUE
    ##_______________Collect plotting bounds info______________
    ##_zoom in scaling with expanding to avoid hexagons outside plot frame__

    if(is(bin1,"erodebin")) {
      tmp1$x <- tmp1$x[bin1@eroded]
      tmp1$y <- tmp1$y[bin1@eroded]
      nxbnds1 <- if(is.null(xbnds1)) range(tmp1$x) else xbnds1
      nybnds1 <- if(is.null(ybnds1)) range(tmp1$y) else ybnds1
      ratiox <- diff(nxbnds1)/diff(bin1@xbnds)
      ratioy <- diff(nybnds1)/diff(bin1@ybnds)
      ratio <- max(ratioy, ratiox)
      nxbnds1 <- mean(nxbnds1) + c(-1,1)*(unzoom * ratio * diff(bin1@xbnds))/2
      nybnds1 <- mean(nybnds1) + c(-1,1)*(unzoom * ratio * diff(bin1@ybnds))/2
    }
    else{
      nxbnds1 <- if(is.null(xbnds1)) bin1@xbnds else xbnds1
      nybnds1 <- if(is.null(ybnds1)) bin1@ybnds else ybnds1
    }

    if(is(bin2,"erodebin")){
      tmp2$x <- tmp2$x[bin2@eroded]
      tmp2$y <- tmp2$y[bin2@eroded]

      nxbnds2 <- if(is.null(xbnds2)) range(tmp2$x) else xbnds2
      nybnds2 <- if(is.null(ybnds2)) range(tmp2$y) else ybnds2
      ratiox <- diff(nxbnds2)/diff(bin2@xbnds)
      ratioy <- diff(nybnds2)/diff(bin2@ybnds)
      ratio <- max(ratioy, ratiox)
      nxbnds2 <- mean(nxbnds2) + c(-1,1)*(unzoom * ratio * diff(bin2@xbnds))/2
      nybnds2 <- mean(nybnds2) + c(-1,1)*(unzoom * ratio * diff(bin2@ybnds))/2
    }
    else{
      nxbnds2 <- if(is.null(xbnds2)) bin2@xbnds else xbnds2
      nybnds2 <- if(is.null(ybnds2)) bin2@ybnds else ybnds2
    }

    xbnds <- range(bin1@xbnds,bin2@xbnds)
    ybnds <- range(bin1@ybnds,bin1@ybnds)
    nxbnds <- range(nxbnds1,nxbnds2)
    nybnds <- range(nybnds1,nybnds2)

    ##__________________ Construct hexagons for bin1 and bin2___________________
    xbins <- bin1@xbins
    shape <- bin1@shape
    hexC1 <- hexcoords(dx = (0.5 * diff(bin1@xbnds))/xbins,
                       dy = (0.5 * diff(bin1@ybnds))/(xbins * shape * sqrt(3)))

    xbins <- bin2@xbins
    shape <- bin1@shape
    hexC2 <- hexcoords(dx = (0.5 * diff(bin2@xbnds))/xbins,
                       dy = (0.5 * diff(bin2@ybnds))/(xbins * shape * sqrt(3)))
    ##__Get unique and ovelapping bin cells from each bin object,ungh__________
    hull1 <- hull2poly(chull(tmp1),tmp1)
    hull2 <- hull2poly(chull(tmp2),tmp2)
    inter.poly<-intersect(hull1,hull2)
    ## Check to see if intersection is NULL
    noint <- TRUE
    inpol1 <- rep(length(inter.poly@pts$x),0)
    inpol2 <- rep(length(inter.poly@pts$x),0)
    if(length(inter.poly@pts)>0){
      ## Find the points in the intersection polygon
      inpol1 <- which.in.poly (tmp1,inter.poly)
      inpol2 <- which.in.poly (tmp2,inter.poly)
      noint <- FALSE
    }


    ##_______________ Full Cell Plotting for Unique Bin1 Cells_________________
    if(any(!inpol1))
        hexpolygon(x = tmp1$x[!inpol1], y= tmp1$y[!inpol1], hexC1,
                   density = density[1], border = border[1], fill = pen[1])

    ##_______________ Full Cell Plotting for Unique Bin2 Cells_________________
    if(any(!inpol2))
        hexpolygon(x = tmp2$x[!inpol2], y= tmp2$y[!inpol2], hexC2,
                   density = density[3], border = border[3], fill = pen[3])

    ##_______________ Full Cell Plotting for Common Cells_____________________
    #if(any(!na1))
    #    hexpolygon(x = tmp1$x[!na1], y= tmp1$y[!na1], hexC,
    #               density = density[2], border = border[2], col = pen[2])
    if(length(inter.poly@pts)>0){
      ipp <- get.pts(inter.poly)
      grid.polygon(ipp[[1]]$x, ipp[[1]]$y,
                   id.lengths=length(ipp[[1]]$x),
                   default.units="native",
                   gp=gpar(col= border[3], fill= pen[2]))
    }

    ##_____________________________Plot Median Cells___________________________
    if(is(bin1,"erodebin") &&is (bin2,"erodebin")) {
        med1 <- which.max(bin1@erode)
        xold <- tmp1$x[med1]
        yold <- tmp1$y[med1]
        med2 <- which.max(bin2@erode)
        xnew <- tmp2$x[med2]
        ynew <- tmp2$y[med2]
        if(abs(xnew - xold) + abs(ynew - yold) > eps) {
            hexpolygon(xold, yold, hexC1,
                       density = density[4], border = border[4], fill = pen[4])
            hexpolygon(xnew, ynew, hexC2,
                       density = density[5], border = border[5], fill = pen[5])
            grid.arrows(c(xold, xnew), c(yold, ynew),
                        length = unit(size,"inches"),
                        default.units="native",
                        gp=gpar(lwd = lwd))
        }
        else {
            hexpolygon(xold, yold, hexC1,
                       density = density[6], border = border[6], fill = pen[6])
        }
    }
    ##________________Clean Up_______________________________________________
    invisible(hvp)
} ## hdiffplot()
