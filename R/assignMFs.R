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
#' assignment <- assignMFs(peakData,p)
#'
#' @export

assignMFs <- function(dat,parameters,verbose = T) {
  options(digits = 10)
  
  if (verbose == T) {
    startTime <- proc.time()
    message(blue('\nMFassign '),red(str_c('v',packageVersion('MFassign') %>% as.character())),' ',date())
    message(rep('_',console_width()))
    params <- parameters %>%
      {capture.output(print(.))} %>%
      {.[-1]} %>%
      {
        .[1] <- yellow(.[1])
        return(.)
      } %>%
      str_c(collapse = '\n')
    message(params)
    message(rep('_',console_width()),'\n')
    message('No. m/z:\t',ncol(dat),'\n')
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
   message(rep('_',console_width()))
   message('\n',green('Complete! '),elapsed,'\n')
 }
 
 return(assignment)
}