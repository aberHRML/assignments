#' assignMFs
#' @description assign molecular formulas to a set of given m/z.
#' @param correlations table containing correlations of m/z to assign molecular formulas
#' @param parameters an S4 object of class AssignmentParamters containing the parameters for molecular formula assignment
#' @param verbose should output be printed to the console
#' @importFrom tibble tibble
#' @importFrom stringr str_split_fixed
#' @importFrom cli console_width
#' @importFrom lubridate seconds_to_period
#' @examples 
#' p <- assignmentParameters('FIE')
#' p@nCores <- 2
#'
#' cors <- correlations[correlations$r > 0.998 | correlations$r < -0.998,]
#'
#' assignment <- assignMFs(cors,p)
#'
#' @export

assignMFs <- function(correlations,parameters,verbose = T) {
  options(digits = 10)
  
  if (verbose == T) {
    startTime <- proc.time()
    cat(blue('\nMFassign'),red(str_c('v',packageVersion('MFassign') %>% as.character())),date(),'\n')
    cat(rep('_',console_width()),'\n',sep = '')
    print(parameters)
    cat(rep('_',console_width()),'\n\n',sep = '')
    cat('Correlations:\t',nrow(correlations),'\n\n')
  }
  
  assignment <- new('Assignment',
                    log = list(date = date(),verbose = verbose),
                    flags = character(),
                    parameters = parameters,
                    correlations = correlations,
                    relationships = tibble(),
                    addIsoAssign = list(),
                    transAssign = list(),
                    assignments  = tibble()
  )
  
 assignment <- assignment %>%
   doAssignment()
 
 if (verbose == T) {
   endTime <- proc.time()
   elapsed <- {endTime - startTime} %>%
     .[3] %>%
     round(1) %>%
     seconds_to_period() %>%
     str_c('[',.,']')
   cat(rep('_',console_width()),'\n',sep = '')
   cat('\n',green('Complete! '),elapsed,'\n\n')
 }
 
 return(assignment)
}