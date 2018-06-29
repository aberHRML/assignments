
assignMethods <- function(method = NULL) {
  methods <- list(
    FIE = FIEassignment,
    `RP-LC` = LCassignment,
    `NP-LC` = LCassignment
  )
  
  if (is.null(method)) {
    method <- methods
  } else {
    method <- methods[[method]]
  }
  return(method)
}