#.RcoclustEnv <- new.env(parent = emptyenv());

.noGenerics <- TRUE

.onLoad=function(libname, pkgname)
  library.dynam("Rcoclust", pkgname, libname)

.onUnload <- function(libpath)
  library.dynam.unload("Rcoclust", libpath)


