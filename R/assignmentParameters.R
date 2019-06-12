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
      p <- new('AssignmentParameters',
          technique = 'FIE',
          maxM = 600,
          maxMFscore = 5,
          ppm = 5,
          limit = 0.001,
          RTwindow = numeric(),
          isotopes = c('13C','18O','13C2'),
          adducts = list(n = c("[M-H]1-", "[M+Cl]1-", "[M+K-2H]1-", 
                               "[M-2H]2-", "[M+Cl37]1-","[2M-H]1-"),
                         p = c('[M+H]1+','[M+K]1+','[M+Na]1+','[M+K41]1+',
                               '[M+NH4]1+','[M+2H]2+','[2M+H]1+')),
          transformations = mzAnnotation::Transformations$`MF Change`,
          nCores = detectCores(),
          clusterType = 'FORK'
      )
    }
    if (technique == 'RP-LC') {
      p <- new('AssignmentParameters',
          technique = 'RP-LC',
          maxM = 600,
          maxMFscore = 5,
          ppm = 5,
          limit = 0.001,
          RTwindow = 1/60,
          isotopes = c('13C','18O','13C2'),
          adducts = list(n = c("[M-H]1-", "[M+Cl]1-", "[M+K-2H]1-", 
                               "[M-2H]2-", "[M+Cl37]1-","[2M-H]1-"),
                         p = c('[M+H]1+','[M+K]1+','[M+Na]1+','[M+K41]1+',
                               '[M+NH4]1+','[M+2H]2+','[2M+H]1+')),
          transformations = mzAnnotation::Transformations$`MF Change`,
          nCores = detectCores(),
          clusterType = 'FORK'
      )
    }
    if (technique == 'NP-LC') {
      p <- new('AssignmentParameters',
               technique = 'NP-LC',
               maxM = 600,
               maxMFscore = 5,
               ppm = 5,
               limit = 0.001,
               RTwindow = 1/60,
               isotopes = c('13C','18O','13C2'),
               adducts = list(n = c('[M-H]1-','[M+Hac-H]1-','[M-2H]2-','[2M-H]1-'),
                              p = c('[M+H]1+','[M+NH4]1+','[M+H-H2O]1+','[M+ACN+H]1+','[M+Na]1+','[M+2H]2+','[2M+H]1+')),
               transformations = mzAnnotation::Transformations$`MF Change`,
               nCores = detectCores(),
               clusterType = 'FORK'
      )
    }
    if (.Platform$OS.type == 'windows') {
      p@clusterType <- 'PSOCK'
    }
    p@nCores <- {detectCores() * 0.75} %>% round()
  }
  return(p)
} 