scattersmooth <- function(x,y, nbin = 200, lambda = 10, ndot = 500, ...)
{

    fillhist <- function(xb, yb, nb)
    {
        H <- matrix(rep(0, prod(nb)), nb[1], nb[2])
        for (i in 1:length(xb))
        {
            H[xb[i],yb[i]] <- H[xb[i],yb[i]] + 1
        }   
        H
    }
 
    # check correct input
    if (length(x) != length(y))
        stop("lengths of x and y do not match")
    if ( (length(x) < 2) | (length(y) < 2) )
        stop("x and y should be vectors")
    if ( !is.numeric(x) | !is.numeric(y) )
        stop("x and y should contain numeric values")
    if ( all(!is.numeric(lambda)) | all(lambda < 0) | (length(lambda) > 2) )
        stop("lambda should be numeric and positive")
    if ( length(lambda) == 1 )
        lambda <- c(lambda, lambda)
    if (all(!is.numeric(nbin)) | all(nbin < 1) | (length(nbin) > 2) )
        stop("nbin should be a strictly positive integer")
    if ( length(nbin) == 1 )
        nbin <- c(nbin, nbin)
    if (!is.numeric(ndot) | (ndot < 0) | (length(ndot) > 1) )
        stop("ndot should be a strictly positive integer")
    ndot <- floor(ndot)
    m <- length(x)

    # Put the x-values into bins
    xmin <- min(x)
    xmax <- max(x)
    dx <- (xmax - xmin) / (nbin[1] - 1)
    xbin <- floor(1 + (x - xmin) / dx)
    xscale <- xmin + (1:nbin[1] - 0.5) * dx

    # Put the y-values into bins
    ymin <- min(y)
    ymax <- max(y)
    dy <- (ymax - ymin) / (nbin[2]-1)
    ybin <- floor(1 + (y - ymin) / dy)
    yscale <- ymin + (1:nbin[2] - 0.5) * dy

    # Create the unsmoothed histogram
    H <- fillhist(xbin, ybin, nbin)

    # Calculate the smoothing matrix
    D1x <- diff(diag(nbin[1]))
    D1y <- diff(diag(nbin[2]))
    D2x <- diff(D1x)
    D2y <- diff(D1y)    
    Qx <- diag(nbin[1]) + lambda[1]^2 * t(D2x) %*% D2x + 2 * lambda[1] * t(D1x) %*% D1x
    Qy <- diag(nbin[2]) + lambda[2]^2 * t(D2y) %*% D2y + 2 * lambda[2] * t(D1y) %*% D1y
    
    # Smooth
    H <- t(solve(Qy, t(solve(Qx, H))))
    
    # Plot coloured image
    image(x = xscale, y = yscale, z = -H, xlab = "", ylab = "", col = heat.colors(100), ...)
    #image(x = xscale, y = yscale, z = -H, xlab = "", ylab = "", col = gray(100), ...)
    
    # Plot selection of dots
    if (ndot > 0) {
      ndot <- min(m, ndot)
      sel <- sort.list(rnorm(m))[1:ndot]
      points(x[sel], y[sel], cex = 0.1) 
    }
}
