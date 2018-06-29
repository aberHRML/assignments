
FIEassignment <- function(element = NULL) {
  methods <- list(
    prepCorrelations = function(assignment){
      assignment %>%
        prepCorrelations()
    },
    relationships = function(assignment){
      assignment %>% 
        relationships()
    },
    adductIsotopeAssignment = function(assignment){
      assignment %>%
        addIsoAssign()
    },
    transformationAssignment = function(assignment){
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
    }
  )
  
  if (!is.null(elements)) {
    return(methods[[element]])
  } else {
    return(methods) 
  }
}