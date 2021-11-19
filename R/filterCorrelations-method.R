
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