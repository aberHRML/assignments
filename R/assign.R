#' Perform molecular formula assignment
#' @rdname assign
#' @description Perform automated molecular formula assignment.
#' @param feature_data a tibble or an object of S4 class `AnalysisData` or `Analysis` containing the feature intensity matrix of m/z for which to assign molecular formulas. See details.
#' @param parameters an S4 object of class `AssignmentParamters` containing the parameters for molecular formula assignment
#' @param verbose should progress output be printed to the console
#' @param type `pre-treated` or `raw` data on which to perform assignment when argument `feature_data` is of S4 class `Analysis`
#' @param ... arguments to pass to the relevant method
#' @details 
#' If argument `feature_data` is specified as a tibble, this should be a feature intensity matrix where 
#' the columns are the `m/z` features to assign and the rows are the individual observations, with the 
#' cells as abundance values. he m/z features provided as column names should be in the form of 
#' <ionisation_mode><m/z>@<retention_time>. Ionisation mode should be given as a prefix n or p for negative 
#' or positive ionisation modes respectively. Feature m/z should be provided to an accuracy of least 5 decimal 
#' places. The retention time portion (@<retention_time>) is only required for LC-MS data and should be provided 
#' in minutes.
#' @importFrom tibble tibble
#' @importFrom stringr str_split_fixed
#' @importFrom cli console_width
#' @importFrom lubridate seconds_to_period
#' @importFrom utils capture.output
#' @examples 
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' assignments <- assignMFs(feature_data,p)
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
    message(blue('\nassignments '),red(str_c('v',packageVersion('assignments') %>% as.character())),' ',date())
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
                    data = feature_data,
                    log = list(verbose = verbose))
  
  assignment <- calcCorrelations(assignment)
  assignment <- calcRelationships(assignment)
  assignment <- addIsoAssign(assignment)
  assignment <- transformationAssign(assignment)
  
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
