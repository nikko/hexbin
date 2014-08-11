panel.hexagons<-function()
{
}

hmatplot<-function()
{
  UseMethod("plot")
}

hmatplot.default <- function(namemat, rlabels, clabels,
                     brd = c(0.5, 0.7, 0.7, 0.5),
                     mai = rep(0.4, 4), unzoom = 1.04, cex = 1.5,
                     border = list(hbox=c(FALSE, FALSE),
                                   hdiff=rep(FALSE,6)),
                     pen = list( hbox = c(2, 3),hdiff = 4:9),
                     density = c(-1, -1, -1, 0, 0, 0),
                     size = 0.06, lwd = 2)
{
    frame() ## <<- MM: always?

    ##___Get initial constants and check arguments_______________________
    nr <- nrow(namemat)
    nc <- ncol(namemat)
    if(missing(rlabels))
        rlabels <- rep("", nr)
    if(missing(clabels))
        clabels <- rep("", nc)
    if(any(mai > brd))
        stop("Border `brd' must be larger than margin `mai' on every side")
    ##_____________________________________________________________
    ##
    ##  Step 1 - Find x and y ranges that maximize the resolution
    ##           and maintain the shape
    tnam <- as.vector(namemat)
    bin1 <- get(tnam[1])
    tmp <- hcell2xy(bin1)
    rx <- range(tmp$x)
    ry <- range(tmp$y)
    for(i in 2:length(tnam)) {
        tmp <- hcell2xy(get(tnam[i]))
        rx <- range(rx, tmp$x)
        ry <- range(ry, tmp$y)
    }
    xbnds <- bin1$xbnds
    ybnds <- bin1$ybnds
    ratiox <- diff(rx)/diff(xbnds)
    ratioy <- diff(ry)/diff(ybnds)
    ratio <- max(ratioy, ratiox)

    rx <- mean(rx) + c(-1,1)*(unzoom * ratio * diff(xbnds))/2
    ry <- mean(ry) + c(-1,1)*(unzoom * ratio * diff(ybnds))/2
    ##_________________________________________________________
    ##  Step 2.  set up a matrix of plots
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    opar<-list()
    opar$global<- par(xpd = TRUE, xaxt = "n", yaxt = "n")
    din <- par("din")
    xsize <- din[1] - brd[2] - brd[4]
    ysize <- din[2] - brd[1] - brd[3]
    nrplot <- 2 * nr - 1
    ncplot <- 2 * nc - 1
    nshape <- (ysize/nrplot)/(xsize/ncplot)
    shape <- bin1$shape
    if(nshape < shape) {
        inc <- (xsize * (1 - nshape/shape))/2
        brd <- brd + c(0, inc, 0, inc)
        xsize <- xsize - 2 * inc
    }
    else {
        inc <- (ysize * (1 - shape/nshape))/2
        brd <- brd + c(inc, 0, inc, 0)
        ysize <- ysize - 2 * inc
    }
    pxsize <- xsize/ncplot
    pysize <- ysize/nrplot
    st <- c((brd[2] - mai[2])/din[1],
            (brd[2] + pxsize + mai[4])/din[1],
            (din[2] - brd[3] - pysize - mai[1])/ din[2],
            (din[2] - brd[3] + mai[3])/din[2])
    pxsize <- pxsize/din[1]
    pysize <- pysize/din[2]
    ##_________________________________________________________
    ## Step 3. produce each boxplot in the structure
    opar$boxpar<-list()
    for(i in seq(length=nr)) {
        for(j in seq(length=nc)) {
            px <- 2 * (j - 1) * pxsize
            py <- 2 * (i - 1) * pysize
            opar$boxpar <- c(opar$boxpar,
                             par(new=TRUE, fig = st + c(px, px, -py, -py),
                                 mai = mai))
            bin <- get(namemat[i, j])
            hboxplot(bin, rx, ry, border=border$hbox, pen = pen$hbox)
            if(i == 1)
                mtext(side = 3, clabels[j], line = 1, cex = cex)
            if(j == 1)
                mtext(side = 2, rlabels[i], line = 1, cex = cex)
        }
    }
    ##_____________________________________________________________
    ## Step 4. produce a difference plot for adjacent pairs in each row
    opar$diffpar.row<-list()
    if(nc > 1) {
        for(i in 1:nr) {
            py <- 2 * (i - 1) * pysize
            for(j in 1:(nc - 1)) {
                px <- (2 * j - 1) * pxsize
                opar$diffpar.row <- c(opar$diffpar.row,
                                      par(new=TRUE,
                                          fig = st + c(px, px, -py, -py)))
                bin1 <- get(namemat[i, j])
                bin2 <- get(namemat[i, j + 1])

                hdiffplot(bin1, bin2, rx, ry, pen = pen$hdiff,
                          border = border$hdiff,
                          size = size, lwd = lwd, xaxt = "n", yaxt = "n",
                          main = "")
            }
        }
    }
    ##_____________________________________________________________
    ## Step 5.  plot difference plot for adjacent pairs in each column
    opar$diffpar.col<-list()
    if(nr > 1) {
        for(j in 1:nc) {
            px <- 2 * (j - 1) * pxsize
            for(i in 1:(nr - 1)) {
                py <- (2 * i - 1) * pysize
                of <- par(new=TRUE, fig = st + c(px, px,  - py,  - py))
                if(!any("fig" == names(opar)))
                  opar$diffpar.col <- c(opar$diffpar.col, of)
                bin1 <- get(namemat[i, j])
                bin2 <- get(namemat[i + 1, j])
                hdiffplot(bin1, bin2, rx, ry, pen = pen$hdiff,
                          border = border$hdiff,
                          size = size, lwd = lwd, xaxt = "n", yaxt = "n",
                          main = "")
            }
        }
    }
    ##_______________________________________________________________
    ## Step 6.  cleanup
    invisible(list(xbnds = rx, ybnds = ry, xsize=xsize, ysize=ysize,
                   plot.par = opar))
    #invisible(list(xbnds = xbnds, ybnds = ybnds, plot.par = opar))
}## hmatplot()

hmatplot.formula <- function(formula, data = parent.frame(), ..., subset, ylab = varnames[response])
{
}
