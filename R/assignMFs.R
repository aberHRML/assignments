#' assignMFs
#' @description assign molecular formulas to a set of given m/z.
#' @param dat tibble containing the peak intensities of m/z for which to assign molecular formulas
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
#' dat <- dplyr::select(peakData,n191.01962,n192.02306,n193.02388,n226.99693,n228.97636,n228.99274,n231.00069,n384.0495,n385.04874)
#'
#' assignment <- assignMFs(dat,p)
#'
#' @export

assignMFs <- function(dat,parameters,verbose = T) {
  options(digits = 10)
  
  if (verbose == T) {
    startTime <- proc.time()
    cat(blue('\nMFassign'),red(str_c('v',packageVersion('MFassign') %>% as.character())),date(),'\n')
    cat(rep('_',console_width()),'\n',sep = '')
    print(parameters)
    cat(rep('_',console_width()),'\n\n',sep = '')
    cat('No. m/z:\t',ncol(dat),'\n\n')
  }
  
  assignment <- new('Assignment',
                    log = list(date = date(),verbose = verbose),
                    flags = character(),
                    parameters = parameters,
                    data = dat,
                    correlations = tibble(),
                    preparedCorrelations = tibble(),
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