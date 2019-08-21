
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
    cat(blue('Filtering correlations '),cli::symbol$continue,'\r',sep = '')
  }
  
  cors <- assignment@correlations %>%
    filterCors(rthresh = assignment@parameters@filter$rthresh,
               n = assignment@parameters@filter$n,
               rIncrement = assignment@parameters@filter$rIncrement,
               nIncrement = assignment@parameters@filter$nIncrement
               )
  
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
    cat(blue('Filtering correlations '),'\t',green(cli::symbol$tick),' ',ncors,' ',elapsed,'\n',sep = '')
  }
  return(assignment)
})