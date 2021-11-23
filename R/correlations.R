
setGeneric('calcCorrelations', function(assignment)
  standardGeneric('calcCorrelations'))

#' @importFrom metabolyseR analysisParameters metabolyse analysisResults keepFeatures analysisData dat sinfo
#' @importFrom magrittr set_rownames
#' @importFrom stats cutree dist hclust

setMethod('calcCorrelations',signature = 'Assignment',function(assignment){
  if (assignment@log$verbose == T) {
    startTime <- proc.time()
    message(blue('Calculating correlations '),cli::symbol$continue,'\r',appendLF = FALSE)
  }
  
  p <- analysisParameters('correlations')
  parameters <- as(assignment,'AssignmentParameters')
  
  p@correlations <- parameters@correlations_parameters
  
  if (str_detect(parameters@technique,'LC')) {
    feat <- tibble(Feature = colnames(assignment@data)) %>%
      mutate(RT = str_split_fixed(Feature,'@',2)[,2] %>%
               as.numeric())
    RTgroups <- feat %>%
      data.frame() %>%
      set_rownames(.$Feature) %>%
      select(-Feature) %>%
      dist() %>%
      hclust() %>%
      cutree(h = parameters@RTwindow) %>%
      {tibble(Feature = names(.),Group = .)}
    
    RTsum <- RTgroups %>%
      group_by(Group) %>%
      summarise(N = n())
    
    RTgroups <- RTgroups %>%
      filter(Group %in% {RTsum %>%
          filter(N > 1) %>%
          .$Group})
    
    cors <- RTgroups %>%
      split(.$Group) %>%
      future_map(~{
        analysisData(assignment@data,tibble(ID = 1:nrow(assignment@data))) %>%
          keepFeatures(features = .x$Feature) %>%
          {metabolyse(dat(.),sinfo(.),p,verbose = FALSE)} %>% 
          analysisResults(element = 'correlations')
      }) %>%
      bind_rows(.id = 'RT Group')
    
  } else {
    cors <- metabolyse(assignment@data,
                       tibble(ID = 1:nrow(assignment@data)),
                       p,
                       verbose = FALSE) %>%
      analysisResults(element = 'correlations') 
  }
  
  assignment@correlations <- cors
  
  if (assignment@log$verbose == T) {
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
                     ID = 1:nrow(.)
              ) %>%
              select(Feature1,Feature2,Mode1:RetentionTime2,log2IntensityRatio,r,ID)
            
            assignment@preparedCorrelations <- correlations
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              message(blue('Preparing correlations '),'\t\t',green(cli::symbol$tick),' ',elapsed)
            }
            
            return(assignment)
          })

setGeneric('filterCorrelations', function(assignment)
  standardGeneric('filterCorrelations'))

setMethod('filterCorrelations',signature = 'Assignment',function(assignment){
  if (assignment@log$verbose == T) {
    startTime <- proc.time()
    message(blue('Filtering correlations '),cli::symbol$continue,'\r',appendLF = FALSE)
  }
  
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
  
  if (assignment@log$verbose == T) {
    endTime <- proc.time()
    elapsed <- {endTime - startTime} %>%
      .[3] %>%
      round(1) %>%
      seconds_to_period() %>%
      str_c('[',.,']')
    ncors <- nrow(assignment@preparedCorrelations) %>%
      str_c('[',.,' correlations',']')
    message(blue('Filtering correlations '),'\t\t',green(cli::symbol$tick),' ',ncors,' ',elapsed)
  }
  return(assignment)
})