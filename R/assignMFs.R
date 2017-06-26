#' @export
#' @importFrom tibble tibble

assignMFs <- function(correlations,parameters) {
  new('Annotation',
      parameters = parameters,
      correlations = correlations,
      relationships = tibble(),
      annotations = tibble(),
      results = tibble()
      )
}