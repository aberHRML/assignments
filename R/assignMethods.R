
assignMethods <- function(method = NULL) {
  methods <- list(
    FIE = function(assignment) {
      assignment <- assignment %>% 
        prepCorrelations() %>%
        relationships() %>% 
        addIsoAssign()
      count <- 0
      while (T) {
        count <- count + 1
        assignment <- suppressWarnings(transformationAssign(assignment))
        if (nrow(assignment@transAssign[[count]]$assigned) == 0) {
          assignment@transAssign <- assignment@transAssign[-count] 
          break()
        }
      }
      return(assignment)
    },
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