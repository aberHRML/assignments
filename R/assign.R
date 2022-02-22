
assignMethods <- function(method) {
  methods <- list(
    FIE = FIEassignment,
    `RP-LC` = LCassignment,
    `NP-LC` = LCassignment
  )
  
  method <- methods[[method]]
  
  return(method)
}

setGeneric('doAssignment',function(assignment)
  standardGeneric('doAssignment')
)

setMethod('doAssignment',signature = 'Assignment',
          function(assignment){
            parameters <- as(assignment,'AssignmentParameters')
            
            assignmentMethod <- assignMethods(parameters@technique)
            
            elements <- names(assignmentMethod())
            elements <- elements[!(elements %in% assignment@flags)]
            
            for(i in elements){
              method <- assignmentMethod(i)
              assignment <- method(assignment)
              assignment@flags <- c(assignment@flags,i)
            }
            
            return(assignment)
          })

#' Assign molecular formulas
#' @rdname assign
#' @description assign molecular formulas to a set of given m/z.
#' @param dat tibble containing the peak intensities of m/z for which to assign molecular formulas
#' @param parameters an S4 object of class AssignmentParamters containing the parameters for molecular formula assignment
#' @param verbose should output be printed to the console
#' @importFrom tibble tibble
#' @importFrom stringr str_split_fixed
#' @importFrom cli console_width
#' @importFrom lubridate seconds_to_period
#' @importFrom utils capture.output
#' @examples 
#' plan(future::sequential)
#' p <- assignmentParameters('FIE')
#'
#' assignment <- assignMFs(peakData,p)
#'
#' @export

assignMFs <- function(dat,parameters,verbose = TRUE) {
  options(digits = 10)
  
  if (verbose == TRUE) {
    startTime <- proc.time()
    message(blue('\nMFassign '),red(str_c('v',packageVersion('MFassign') %>% as.character())),' ',date())
    message(rep('_',console_width()))
    params <- parameters %>%
      {capture.output(print(.))} %>%
      {.[-1]} %>%
      {
        .[1] <- yellow(.[1])
        .
      } %>%
      str_c(collapse = '\n')
    message(params)
    message(rep('_',console_width()),'\n')
    message('No. m/z:\t',ncol(dat),'\n')
  }
  
  assignment <- new('Assignment',
                    parameters,
                    data = dat)
  assignment@log$verbose <- verbose
  
  assignment <- assignment %>%
    doAssignment()
  
  if (verbose == TRUE) {
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

#' @rdname assign
#' @export

setGeneric('continueAssignment',function(assignment)
  standardGeneric('continueAssignment')
)


setMethod('continueAssignment',signature = 'Assignment',
          function(assignment){
            doAssignment(assignment)
          })