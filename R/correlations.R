
setGeneric('calcCorrelations', function(assignment)
  standardGeneric('calcCorrelations'))

#' @importFrom metabolyseR analysisParameters metabolyse analysisResults

setMethod('calcCorrelations',signature = 'Assignment',function(assignment){
  if (assignment@log$verbose == T) {
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
                     log2IntensityRatio,r,ID)
            
            assignment@correlations <- correlations
            
            return(assignment)
          })
