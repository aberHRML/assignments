#' assignMFs
#' @description assign molecular formulas to a set of given m/z.
#' @param correlations table containing correlations of m/z to assign molecular formulas
#' @param parameters an S4 object of class AssignmentParamters containing the parameters for molecular formula assignment
#' @importFrom tibble tibble
#' @importFrom stringr str_split_fixed
#' @examples 
#' \dontrun{
#' res <- assignMFs(correlations,assignmentParameters('FIE'))
#' }
#' @export

assignMFs <- function(correlations,parameters) {
  options(digits = 10)
  
  assignment <- new('Assignment',
                    flags = character(),
                    parameters = parameters,
                    correlations = correlations,
                    relationships = tibble(),
                    addIsoAssign = list(),
                    transAssign = list(),
                    assignments  = tibble()
  )
  
 assignment %>%
   doAssignment()
}