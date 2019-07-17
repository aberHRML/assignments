
FIEassignment <- function(element = NULL) {
  methods <- list(
    `calculate correlations` = function(assignment){
      assignment %>%
        calcCorrelations()
    },
    `filter correlations` = function(assignment){
      assignment %>%
        filterCorrelations()
    },
    `prepare correlations` = function(assignment){
      assignment %>%
        prepCorrelations()
    },
    relationships = function(assignment){
      assignment %>% 
        relationships()
    },
    `adduct and isotope assignment` = function(assignment){
      assignment %>%
        addIsoAssign()
    },
    `transformation assignment` = function(assignment){
      count <- 0
      while (T) {
        count <- count + 1
        assignment <- suppressWarnings(transformationAssign(assignment))
        if (length(assignment@transAssign[[count]]) == 0) {
          assignment@transAssign <- assignment@transAssign[-count] 
          break()
        }
        if (nrow(assignment@transAssign[[count]]$assigned) == 0) {
          assignment@transAssign <- assignment@transAssign[-count]  
          break()
        }
         
      }
      return(assignment)
    }
  )
  
  if (!is.null(element)) {
    return(methods[[element]])
  } else {
    return(methods) 
  }
}