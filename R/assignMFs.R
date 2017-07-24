#' @importFrom tibble tibble
#' @export

assignMFs <- function(correlations,parameters) {
  annotation <- new('Annotation',
      parameters = parameters,
      correlations = correlations,
      relationships = tibble(),
      addIsoAssign = list(),
      transAssign = list(),
      assignments  = tibble()
      )
  method <- assignMethods(annotation@parameters@technique)
  annotation <- method(annotation)
  return(annotation)
}