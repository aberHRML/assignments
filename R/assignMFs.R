#' @importFrom tibble tibble
#' @export

assignMFs <- function(correlations,parameters) {
  assignment <- new('Assignment',
      parameters = parameters,
      correlations = correlations,
      relationships = tibble(),
      addIsoAssign = list(),
      transAssign = list(),
      assignments  = tibble()
      )
  method <- assignMethods(assignment@parameters@technique)
  assignment <- method(assignment)
  return(assignment)
}