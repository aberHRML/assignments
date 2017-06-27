#' @export
#' @importFrom tibble tibble

assignMFs <- function(correlations,parameters) {
  new('Annotation',
      parameters = parameters,
      correlations = correlations,
      relationships = tibble(),
      addIsoAssign = list(),
      transAssign = list(),
      assignments  = tibble()
      )
}