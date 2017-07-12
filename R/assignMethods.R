
assignMethods <- function(method = NULL) {
  methods <- list(
    FIE = function(annotation) {
      annotation <- annotation %>% 
        relationships() %>% 
        addIsoAssign()
      return(annotation)
    }
  )
  
  if (is.null(method)) {
    method <- methods
  } else {
    method <- methods[[method]]
  }
  return(method)
}