grid.hexagon <- function (x = 0.5, y = 0.5, r = 0.5, s = 1,
                          default.units = "npc",
                          name = NULL, gp = gpar(), draw = TRUE, vp = NULL) 
{
    hg <- hexagonGrob(x = x, y = y, r = r, s=s, default.units = default.units, 
        name = name, gp = gp, vp = vp)
    if (draw) 
        grid.draw(hg)
    invisible(hg)
}

hexagonGrob <- 
function (x = 0.5, y = 0.5, r = 0.5, s = 1, default.units = "npc", name = NULL, 
          gp = gpar(), vp = NULL) 
{
    if (!is.unit(x)) 
        x <- unit(x, default.units)
    if (!is.unit(y)) 
        y <- unit(y, default.units)
    if (!is.unit(r)) 
        r <- unit(r, default.units)
    if (!is.numeric(s))
      s <- as.numeric(s)
    grob(x = x, y = y, r = r, s=s, name = name, gp = gp, vp = vp, 
        cl = "hexagon")
 
}
######################################
# HEXAGON primitive
######################################
validDetails.hexagon <- function(x) {
  if (!is.unit(x$x) ||
      !is.unit(x$y) ||
      !is.unit(x$r))
    stop("'x', 'y', and 'r' must be units")
  if (!is.numeric(s))
    stop("shape parameter 's', must be numeric")
  x
}



drawDetails.hexagon <- function(x, recording=TRUE) {
  
  
  grid.Call.graphics("L_polygon", h$x, h$y,
                       list(as.integer(1:length(x$x))))
}

xDetails.polygon <- function(x, theta) {
    bounds <- grid.Call("L_locnBounds", x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[1], "inches")
}

yDetails.polygon <- function(x, theta) {
    bounds <- grid.Call("L_locnBounds", x$x, x$y, theta)
    if (is.null(bounds))
        unit(0.5, "npc")
    else
        unit(bounds[2], "inches")
}

widthDetails.polygon <- function(x) {
  bounds <- grid.Call("L_locnBounds", x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[3], "inches")
}

heightDetails.polygon <- function(x) {
  bounds <- grid.Call("L_locnBounds", x$x, x$y, 0)
  if (is.null(bounds))
    unit(0, "inches")
  else
    unit(bounds[4], "inches")
}

polygonGrob <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                        id=NULL, id.lengths=NULL,
                        default.units="npc",
                        name=NULL, gp=gpar(), vp=NULL) {
  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)
  grob(x=x, y=y, id=id,
       id.lengths=id.lengths,
       name=name, gp=gp, vp=vp, cl="polygon")
}

grid.polygon <- function(x=c(0, 0.5, 1, 0.5), y=c(0.5, 1, 0.5, 0),
                         id=NULL, id.lengths=NULL,
                         default.units="npc",
                         name=NULL, gp=gpar(), draw=TRUE, vp=NULL) {
  pg <- polygonGrob(x=x, y=y, id=id, id.lengths=id.lengths,
                    default.units=default.units,
                    name=name, gp=gp, vp=vp)
  if (draw)
    grid.draw(pg)
  invisible(pg)
}
