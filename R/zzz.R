.onLoad <- function(libname, pkgname) {
  if (options()$digits < 10) options(digits = 10)
  
  invisible()
}