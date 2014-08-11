#Lib load function
.First.lib <- function(libname, pkgname, where) {
  require(grid)
  require(methods)
  require(colorspace)
  library.dynam("hexbin", pkgname, libname)
  #where <- match(paste("package:", pkgname, sep=""), search())
  #.initClasses(where)

}
