#' annotationParameters
#' @description Return assignment parameters for a specified technique.
#' @param technique technique to use for assignment. \code{NULL} prints available techniques
#' @importFrom parallel detectCores
#' @importFrom methods new
#' @export

assignmentParameters <- function(technique = NULL){
  availTechniques <- c('FIE','RP-LC','NP-LC')
  if (is.null(technique)) {
    cat('\nAvailable Techniques:',str_c('\n\t\t\t',str_c(availTechniques,collapse = '\n\t\t\t'),'\n'))
    p <- NULL
  } else {
    
    if (technique == 'FIE') {
      p <- new('AssignmentParameters')
    }
    if (technique == 'RP-LC') {
      p <- new('AssignmentParameters',
          technique = 'RP-LC',
          RTwindow = 1/60
      )
    }
    if (technique == 'NP-LC') {
      p <- new('AssignmentParameters',
               technique = 'NP-LC',
               RTwindow = 1/60
      )
    }
  }
  return(p)
} 