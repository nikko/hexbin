hgridcent<- function(xbins,xbnds,ybnds,shape)
{
    jmax <- floor(xbins + 1.5001)
    imax <- 2 * floor((xbins * hb@shape)/sqrt(3) + 1.5001)
    dimen <- c(imax, jmax)
    c3 <- diff(xbnds)/xbins
    c4 <- (diff(ybnds) * sqrt(3))/(2 * shape * xbins)
    jmax <- dimen[2]
    cell <- 1:(dimen[1]*dimen[2])-1
    i <- cell %/% jmax
    j <- cell %% jmax
    y <- c4 * i + ybnds[1]
    x <- c3 * ifelse(i %% 2 == 0, j, j + 0.5) + xbnds[1]
    list(x = x, y = y, dimen=dimen)
}
interp.2d <- function(x,y,z,x0,y0)
{
	# Bilinear interpolation of function (x,y,z) onto (x0,y0).
	# Taken from Numerical Recipies (second edition) section 3.6.
	# Slow because of loop.  Would be better in C or Fortran.
	# Called by hdr.info.2d

	n <- length(x)
	n0 <- length(x0)
	z0 <- numeric(length=n0)
        
	for(i in 1:n0)
	{
		d <- x - x0[i]
		j <- (1:(n-1))[ d[1:(n-1)]*d[2:n] < 0]
		d <- y - y0[i]
		k <- (1:(n-1))[ d[1:(n-1)]*d[2:n] < 0]
		v <- (x0[i] - x[j])/(x[j+1]-x[j])
		u <- (y0[i] - y[k])/(y[k+1]-y[k])
		z0[i] <- (1-v)*(1-u)*z[j,k] +
                          v*(1-u)*z[j+1,k] +
                          v*u*z[j+1,k+1] +
                          (1-v)*u*z[j,k+1]
	}
	return(z0)
}

interp.hexbin<-function(hbin,method=c("bilinear","cubic"),xbins,
                        xrange=NULL,yrange=NULL)
{
  if(is.null(xrange)) xrange <- hb@xbnds
  if(is.null(yrange)) yrange <- hb@ybnds
  method <- match.arg(method)
  xy <- hgridcent(xbins,xrange,yrange,hb@shape)
  if(method=="bilinear")
    count <- interp.2d(hbin@xcm,hbin@ycm,hbin@count,xy$x,xy$y)
  if(method=="cubic"){
    if(require(akima)){
      count <- interpp(hbin@xcm,hbin@ycm,hbin@count,xy$x,xy$y,ncp=6)$z
      keep <- count > .95 & !is.na(count)
    }
  }
  cell <- 1:(jmax*imax)

  new("hexbin",
      cell = cell[keep], count = count[keep],
      xcm = xy$x[keep], ycm = xy$y[keep], xbins = xbins,
      shape = hb@shape, xbnds = xrange , ybnds = yrange,
      dimen = xy$dimen, n = hb@n, ncells = as.integer(sum(keep)),
      call = hb@call, xlab = hb@xlab, ylab = hb@ylab, cID = NULL)
}
