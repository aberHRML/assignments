#' Molecular formula assignment methods
#' @rdname assignment-methods
#' @description These methods provide the access to performing the individual steps of the molecular 
#' formula assignment approach. See Details for more information of when it is best to use these 
#' instead of `assignMFs()`. 
#' @param assignment an object of S4 class `Assignment`
#' @details 
#' In circumstances where the molecular formula assignment approach has high memory requirements, 
#' such as where there are many correlations (> 2 million) or many high *m/z* (>700), it may be 
#' preferable to perform the assignment steps separately as opposed to using `assignMFs()`. This 
#' can reduce the memory overheads required to successfully assign molecular formulas to the data 
#' and also enable the possibility of objects to be saved and/or unloaded between the assignment 
#' steps where needed.
#' @return An object of S4 class `Assignment` containing molecular formula assignments.
#' @examples 
#' \dontrun{
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' mf_assignments <- assignment(feature_data,p)
#' 
#' mf_assignments <- mf_assignments %>% 
#'    calcCorrelations() %>% 
#'    calcRelationships() %>% 
#'    addIsoAssign() %>% 
#'    transformationAssign()
#' }
#' @export

setGeneric('calcCorrelations', function(assignment)
  standardGeneric('calcCorrelations'))

#' @rdname assignment-methods
#' @importFrom metabolyseR analysisParameters metabolyse analysisResults

setMethod('calcCorrelations',signature = 'Assignment',function(assignment){
  if (assignment@log$verbose == TRUE) {
    startTime <- proc.time()
    message(blue('Calculating correlations '),cli::symbol$continue,'\r',appendLF = FALSE)
  }
  
  p <- analysisParameters('correlations')
  parameters <- as(assignment,'AssignmentParameters')
  
  p@correlations <- parameters@correlations_parameters[c('method',
                                                         'pAdjustMethod',
                                                         'corPvalue',
                                                         'minCoef',
                                                         'maxCor')]
  
  assignment@correlations <- metabolyse(assignment@data,
                                        tibble(ID = 1:nrow(assignment@data)),
                                        p,
                                        verbose = FALSE) %>%
    analysisResults(element = 'correlations') 
  
  assignment <- assignment %>% 
    prepCorrelations()
  
  if (assignment@log$verbose == TRUE) {
    endTime <- proc.time()
    elapsed <- {endTime - startTime} %>%
      .[3] %>%
      round(1) %>%
      seconds_to_period() %>%
      str_c('[',.,']')
    ncors <- nrow(assignment@correlations) %>%
      str_c('[',.,' correlations',']')
    message(blue('Calculating correlations '),'\t',green(cli::symbol$tick),' ',ncors,' ',elapsed)
  }
  
  return(assignment)
})

setGeneric('prepCorrelations', function(assignment)
  standardGeneric('prepCorrelations'))

setMethod('prepCorrelations',signature = 'Assignment',
          function(assignment){
            
            correlations <- assignment@correlations
            
            correlations <- correlations %>%
              mutate(Mode1 = str_split_fixed(Feature1,'@',2) %>% 
                       .[,1] %>%
                       str_sub(1,1),
                     Mode2 = str_split_fixed(Feature2,'@',2) %>% 
                       .[,1] %>%
                       str_sub(1,1),
                     `m/z1` = str_split_fixed(Feature1,'@',2) %>% 
                       .[,1] %>% 
                       str_replace_all('[:alpha:]','') %>% 
                       as.numeric(),
                     `m/z2` = str_split_fixed(Feature2,'@',2) %>% 
                       .[,1] %>% 
                       str_replace_all('[:alpha:]','') %>% 
                       as.numeric(),
                     RetentionTime1 = str_split_fixed(Feature1,'@',2) %>% 
                       .[,2] %>%
                       as.numeric(),
                     RetentionTime2 = str_split_fixed(Feature2,'@',2) %>% 
                       .[,2] %>%
                       as.numeric(),
                     RetentionTimeDiff = abs(RetentionTime1 - RetentionTime2),
                     ID = 1:nrow(.)
              ) %>%
              select(Feature1,Feature2,
                     Mode1:RetentionTimeDiff,
                     log2IntensityRatio,coefficient,ID)
            
            assignment@correlations <- correlations
            
            return(assignment)
          })
