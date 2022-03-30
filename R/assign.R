#' Assign molecular formulas
#' @rdname assign
#' @description assign molecular formulas to a set of given m/z.
#' @param feature_data a tibble or an object of S4 class `AnalysisData` or `Analysis` containing the feature intensity matrix of m/z for which to assign molecular formulas. See details.
#' @param parameters an S4 object of class `AssignmentParamters` containing the parameters for molecular formula assignment
#' @param verbose should progress output be printed to the console
#' @param type `raw` or `pre-treated` data on which to perform assignment when argument `feature_data` is of class `Analysis`
#' @param ... arguments to pass to the relevant method
#' @details 
#' If argument `feature_data` is specified as a tibble, this should be a feature intensity matrix where the columns are the `m/z` features to assign and the rows are the individual observations, with the cells as abundance values.
#' @importFrom tibble tibble
#' @importFrom stringr str_split_fixed
#' @importFrom cli console_width
#' @importFrom lubridate seconds_to_period
#' @importFrom utils capture.output
#' @examples 
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' assignment <- assignMFs(feature_data,p)
#'
#' @export

setGeneric('assignMFs',function(feature_data,
                                parameters = assignmentParameters('FIE-HRMS'),
                                verbose = TRUE,
                                ...)
  standardGeneric('assignMFs')
)

#' @rdname assign

setMethod('assignMFs',signature = 'tbl_df',
          function(feature_data,
                   parameters = assignmentParameters('FIE-HRMS'),
                   verbose = TRUE) {
  
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
    message('No. m/z:\t',ncol(feature_data),'\n')
  }
  
  assignment <- new('Assignment',
                    parameters,
                    data = feature_data)
  assignment@log$verbose <- verbose
  
  assignment <- assignment %>% 
    calcCorrelations() %>% 
    relationships() %>% 
    addIsoAssign()
  
  count <- 0
  while (TRUE) {
    count <- count + 1
    assignment <- suppressWarnings(transformationAssign(assignment))
    if (length(assignment@transAssign[[count]]) == 0) {
      assignment@transAssign <- assignment@transAssign[-count] 
      break()
    }
    if (nrow(assignment@transAssign[[count]]$assigned) == 0) {
      assignment@transAssign <- assignment@transAssign[-count]  
      break()
    }
    
  }
  
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
})

#' @rdname assign
#' @importFrom metabolyseR dat

setMethod('assignMFs',signature = 'AnalysisData',
          function(feature_data,
                   parameters = assignmentParameters('FIE'),
                   verbose = TRUE){
            feature_data %>% 
              dat() %>% 
              assignMFs(parameters = parameters,
                        verbose = verbose)
          })

#' @rdname assign
#' @importFrom metabolyseR raw preTreated

setMethod('assignMFs',signature = 'Analysis',
          function(feature_data,
                   parameters = assignmentParameters('FIE'),
                   verbose = TRUE,
                   type = c('pre-treated','raw')){
            
            type <- match.arg(type,
                              choices = c('pre-treated','raw'))
            
            if (type == 'raw') feature_data <- raw(feature_data)
            if (type == 'pre-treated') feature_data <- preTreated(feature_data)
            
            feature_data %>% 
              assignMFs(parameters = parameters,
                        verbose = verbose)
          })
