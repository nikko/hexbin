
hexInfill <- function(hb){
    ##function (xbins, xbnds, ybnds, shape, edge.add = 0)
    xbins <- hb@xbins
    shape <- hb@shape
    xbnds <- hb@xbnds
    ybnds <- hb@ybnds
    dimen <- hb@dimen

    #jmax <- floor(xbins + 1.5001)
    #c1 <- 2 * floor((xbins * shape)/sqrt(3) + 1.5001)
    #imax <- (jmax * c1 - 1)/jmax + 1
    #dimen <- c(imax, jmax)
    c3 <- diff(xbnds)/xbins
    c4 <- (diff(ybnds) * sqrt(3))/(2 * shape * xbins)

    jmax <- dimen[2]
    cell <- 1:(dimen[1] * dimen[2])
    i <- cell%/%jmax
    j <- cell%%jmax
    y <- c4 * i + ybnds[1]
    x <- c3 * ifelse(i%%2 == 0, j, j + 0.5) + xbnds[1]
    list(x = x, y = y, dimen = dimen, dx = c3, dy = c4)
}



