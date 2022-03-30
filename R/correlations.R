
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
  
  p@correlations <- parameters@correlations_parameters
  
  assignment@correlations <- metabolyse(assignment@data,
                     tibble(ID = 1:nrow(assignment@data)),
                     p,
                     verbose = FALSE) %>%
    analysisResults(element = 'correlations') 
  
  assignment <- assignment %>% 
    filterCorrelations() %>% 
    prepCorrelations()
  
  if (assignment@log$verbose == TRUE) {
    endTime <- proc.time()
    elapsed <- {endTime - startTime} %>%
      .[3] %>%
      round(1) %>%
      seconds_to_period() %>%
      str_c('[',.,']')
    ncors <- nrow(assignment@preparedCorrelations) %>%
      str_c('[',.,' correlations',']')
    message(blue('Calculating correlations '),'\t',green(cli::symbol$tick),' ',ncors,' ',elapsed)
  }
  
  return(assignment)
})


filterCors <- function(correlations, rthresh = 0.7, n = 100000, rIncrement = 0.01, nIncrement = 20000){
  filCors <- function(cors,rthresh,n){
    while (nrow(cors) > n) {
      cors <- correlations %>%
        filter(r > rthresh | r < -rthresh)
      rthresh <- rthresh + rIncrement
    }
    return(cors)
  }
  
  while (TRUE) {
    cors <- filCors(correlations,rthresh,n)
    if (nrow(cors) > 0) {
      break()
    } else {
      n <- n + nIncrement
    }
  }
  return(cors)
}

setGeneric('prepCorrelations', function(assignment)
  standardGeneric('prepCorrelations'))

setMethod('prepCorrelations',signature = 'Assignment',
          function(assignment){
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              message(blue('Preparing correlations '),cli::symbol$continue,'\r',appendLF = FALSE)
            }
            
            correlations <- assignment@preparedCorrelations
            
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
            
            assignment@preparedCorrelations <- correlations
            
            return(assignment)
          })

setGeneric('filterCorrelations', function(assignment)
  standardGeneric('filterCorrelations'))

setMethod('filterCorrelations',signature = 'Assignment',function(assignment){

  parameters <- as(assignment,'AssignmentParameters')
  
  if (str_detect(parameters@technique,'LC')) {
    cors <- assignment@correlations %>%
      filter(r < -(parameters@filter$rthresh) | r > parameters@filter$rthresh)
  } else {
    cors <- assignment@correlations %>%
      filterCors(rthresh = parameters@filter$rthresh,
                 n = parameters@filter$n,
                 rIncrement = parameters@filter$rIncrement,
                 nIncrement = parameters@filter$nIncrement
      ) 
  }
  
  assignment@preparedCorrelations <- cors
  
  return(assignment)
})