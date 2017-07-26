#' annotationParameters
#' @description Return assignment parameters for a specified technique.
#' @param technique technique to use for assignment. \code{NULL} prints available techniques
#' @importFrom parallel detectCores
#' @importFrom methods new
#' @export

assignmentParameters <- function(technique = NULL){
  if (is.null(technique)) {
    cat('\nAvailable Techniques:','\t FIE')
  }
  
  if (technique == 'FIE') {
    new('AssignmentParameters',
        technique = 'FIE',
        maxM = 400,
        maxMFscore = 5,
        ppm = 6,
        limit = 0.001,
        isotopes = c('C13','O18','2C13'),
        adducts = list(n = c("[M-H]1-", "[M+Cl]1-", "[M+K-2H]1-", 
                             "[M-2H]2-", "[M+Cl37]1-","[M+K41-2H]1-",
                             "[M+Na-2H]1-","[2M-H]1-","[3M-H]1-",
                             "[M-3H]3-"),
                       p = c('[M+K]1+','[M+H]1+','[M+Na]1+','[M+K41]1+',
                              '[M+NH4]1+','[M+2H]2+','[2M+H]1+','[M+3H]3+')),
        nCores = detectCores(),
        clusterType = 'FORK'
        )
  }
} 