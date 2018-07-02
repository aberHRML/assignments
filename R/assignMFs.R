#' assignMFs
#' @description assign molecular formulas to a set of given m/z.
#' @param correlations table containing correlations of m/z to assign molecular formulas
#' @param parameters an S4 object of class AssignmentParamters containing the parameters for molecular formula assignment
#' @param verbose should output be printed to the console
#' @importFrom tibble tibble
#' @importFrom stringr str_split_fixed
#' @importFrom cli console_width
#' @examples 
#' \dontrun{
#' res <- assignMFs(correlations,assignmentParameters('FIE'))
#' }
#' @export

assignMFs <- function(correlations,parameters,verbose = T) {
  options(digits = 10)
  
  if (verbose == T) {
    cat(blue('\nMFassign'),red(str_c('v',packageVersion('MFassign') %>% as.character())),date(),'\n')
    cat(rep('_',console_width()),'\n')
  }
  
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