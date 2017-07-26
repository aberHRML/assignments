#' assignMFs
#' @description assign molecular formulas to a set of given m/z.
#' @param correlations table containing correlations of m/z to assign molecular formulas
#' @param parameters an S4 object of class AssignmentParamters containing the parameters for molecular formula assignment
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