hcell2xyInt <- function(xbins, xbnds, ybnds, shape)
{  
  jmax <- floor(xbins + 1.5001)
  imax <- 2 * floor((xbins *shape)/sqrt(3) + 1.5001)
  dimen <- c(imax, jmax)
  cell <- 1:(dimen[1]*dimen[2])-1
  i <- cell %/% jmax
  j <- cell %% jmax
  list(i=i+1, j=j+1)
}

xyInt2hcell <- function(){

}

hgridcent <- function(xbins, xbnds, ybnds, shape)
{
  jmax <- floor(xbins + 1.5001)
  c1 <- 2 * floor((xbins *shape)/sqrt(3) + 1.5001)
  imax <- (jmax*c1 -1)/jmax + 1
  dimen <- c(imax, jmax)
  c3 <- diff(xbnds)/xbins
  c4 <- (diff(ybnds) * sqrt(3))/(2 * shape * xbins)
  jmax <- dimen[2]
  cell <- 1:(dimen[1]*dimen[2])
  i <- cell %/% jmax
  j <- cell %% jmax
  y <- c4 * i + ybnds[1]
  x <- c3 * ifelse(i %% 2 == 0, j, j + 0.5) + xbnds[1]
  list(x = x, y = y, dimen = dimen, dx=c3, dy=c4)
}

hgridcent2 <- function(dimen,dx,dy,xbnds,ybnds)
{
  jmax <- dimen[2]
  cell <- 1:(dimen[1]*dimen[2]+2)-1
  i <- cell %/% jmax
  j <- cell %% jmax
  y <- dy * i + ybnds[1]
  x <- dx* ifelse(i %% 2 == 0, j, j + 0.5) + xbnds[1]
  list(x = x, y = y, dimen = dimen, dx=dx, dy=dy)
}

hexGraphPaper <- function(hb, xbnds=NULL, ybnds=NULL, xbins=30, shape=1,
                          add=TRUE, fill.edges=1, fill=0, border=1)
{
  if(missing(hb) && (is.null(xbnds) || is.null(ybnds)))
	  stop("Need a hexbin object or boundaries to make lattice")
  if(!missing(hb)){
    xbins <- hb@xbins
    shape <- hb@shape
    xbnds <- if(is.null(xbnds)) hb@xbnds else xbnds
    ybnds <- if(is.null(ybnds)) hb@ybnds else ybnds
    dimen <- hb@dimen
  }
  if(fill.edges<=0)
    xy <- hgridcent(xbins, xbnds, ybnds, shape)
  if((fill.edges>0) || add) {
    sx <- xbins/diff(xbnds)
    sy <- (xbins * shape)/diff(ybnds)
    inner <- 0.5
    outer <- (2 * inner)/sqrt(3)
    dx <- inner/sx
    dy <- outer/(2 * sy)
    if(fill.edges>0) {
      naddx <- fill.edges*1.5*(2*dx)
      naddy <- fill.edges*(2*dy)
      #xy <- hgridcent(xbins+2*fill.edges, xbnds+c(-naddx,naddx),
                      #ybnds+c(-naddy,naddy), shape)
      xy <- hgridcent2(xy$dimen+2, 2*dx, 2*dy, xbnds+c(0,naddx),
                      ybnds+c(0,naddy))
    } 
    if(add) {
      hexC <- hexcoords(dx, dy, sep=NULL)
      hexpolygon (xy$x, xy$y, hexC, dx, dy,
                  fill = fill, border = border, hUnit = "native")
    }
  }
  invisible(xy)
}

hexTapply <- function(hbin,dat,FUN=sum,...)
{
  if(is.null(hbin@cID))
    stop("Must have cell ID's to do this operation \n
          please re-bin data using IDs = TRUE")
  if((length(dat)> 0) & (length(dat) != length(hbin@cID)))
    stop("Length of IDs does not match the length of the data")
  tapply(dat,hbin@cID,FUN,...)
}

hexWapply <- function(hbin,nh.order=1,data=NULL,FUN=sum,...)
{
  #first convert to xy coordinates
  if(!is.null(data) && is.null(hbin@cID))
    stop("Must have cell ID's to do this operation \n
          please re-bin data using IDs = TRUE")
  if(!is.null(data)&&((length(dat)> 0) & (length(dat) != length(hbin@cID))))
    stop("Length of IDs does not match the length of the data")
  xy <- hcell2xyInt(hbin@xbins, hbin@xbnds, hbin@ybnds, hbin@shape){
}

running <- 
function (X, Y = NULL, fun = mean, width = min(length(X), 20),
    allow.fewer = FALSE, pad = FALSE, align = c("right", "center",
        "left"), simplify = TRUE, by, ...)
{
    align <- match.arg(align)
    n <- length(X)
    if (align == "left") {
        from <- 1:n
        to <- pmin((1:n) + width - 1, n)
    }
    else if (align == "right") {
        from <- pmax((1:n) - width + 1, 1)
        to <- 1:n
    }
    else {
        from <- pmax((2 - width):n, 1)
        to <- pmin(1:(n + width - 1), n)
        if (!odd(width))
            stop("width must be odd for center alignment")
    }
    elements <- apply(cbind(from, to), 1, function(x) seq(x[1],
        x[2]))
    if (is.matrix(elements))
        elements <- as.data.frame(elements)
    names(elements) <- paste(from, to, sep = ":")
    if (!allow.fewer) {
        len <- sapply(elements, length)
        skip <- (len < width)
    }
    else {
        skip <- 0
    }
    run.elements <- elements[!skip]
    if (!invalid(by))
        run.elements <- run.elements[seq(from = 1, to = length(run.elements),
            by = by)]
    if (is.null(Y)) {
        funct <- function(which, what, fun, ...) fun(what[which],
            ...)
        if (simplify)
            Xvar <- sapply(run.elements, funct, what = X, fun = fun,
                ...)
        else Xvar <- lapply(run.elements, funct, what = X, fun = fun,
            ...)
    }
    else {
        funct <- function(which, XX, YY, fun, ...) fun(XX[which],
            YY[which], ...)
        if (simplify)
            Xvar <- sapply(run.elements, funct, XX = X, YY = Y,
                fun = fun, ...)
        else Xvar <- lapply(run.elements, funct, XX = X, YY = Y,
            fun = fun, ...)
    }
    if (allow.fewer || !pad)
        return(Xvar)
    if (simplify)
        if (is.matrix(Xvar)) {
            wholemat <- matrix(new(class(Xvar[1, 1]), NA), ncol = length(to),
                nrow = nrow(Xvar))
            colnames(wholemat) <- paste(from, to, sep = ":")
            wholemat[, -skip] <- Xvar
            Xvar <- wholemat
        }
        else {
            wholelist <- rep(new(class(Xvar[1]), NA), length(from))
            names(wholelist) <- names(elements)
            wholelist[names(Xvar)] <- Xvar
            Xvar <- wholelist
        }
    return(Xvar)
}

hexGapply <- function(hbin,Margin,FUN=sum,...)
{
  stopifnot(is.numeric(Margin) || ( Margin<=3 & Margin>=1 ),)
  if(Margin==1){
    
  }
}

hexTribble <- function(hbin,cc=TRUE,n.tol=2)
{

}

## Finish this function for the next round of hdiffplot
intersect.hexbins <- function(bin1,bin2)
{
  ## For hexbin objects which do not have the same range
  ## and number of bins this function returns the bins
  ## in 1 only, the bins in 2 only, the intersection polygon
  ## and 2 fringes. Fringe 1 is the pieces of bin1 from the hexagons
  ## that partially intersect bin2 and fringe 2 is the same but the
  ## other way around.
  ## 

}

optShape <- function(vp, height=NULL, width=NULL, mar=NULL)
{
  if(missing(vp) && (is.null(height) || is.null(width)))
    stop("Need a viewport object or height and width of the plotting region.") 
  if(!missing(vp)){
    if("hexVP"%in%class(vp)){
      height <- vp@plt[2]
      width <- vp@plt[1]
    }
    else if("viewport"%in%class(vp)){
      #height <- convertHeight(unit(1,"npc"),"inches")
      #width <- convertWidth (unit(1,"npc"),"inches")
      height <- convertUnit(vp$height,"inches")
      width <- convertUnit(vp$width,"inches")
    }
    else stop("need valid viewport or hexViewport")
  }
  if(!is.null(mar)){
    height <- height-mar[1]-mar[3]
    width <- width -mar[2]-mar[4]
  }
  shape <- height/width
  shape
}  
