
assignMethods <- function(method = NULL) {
  methods <- list(
    FIE = function(annotation) {
      annotation <- annotation %>% 
        relationships() %>% 
        addIsoAssign()
      count <- 0
      while (T) {
        count <- count + 1
        annotation <- suppressWarnings(transformationAssign(annotation))
        if (nrow(annotation@transAssign[[count]]$assigned) == 0) {
          annotation@transAssign <- annotation@transAssign[-count] 
          break()
        }
      }
      
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