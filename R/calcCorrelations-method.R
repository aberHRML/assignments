#' @importFrom metabolyseR analysisParameters metabolyse correlationResults

setMethod('calcCorrelations',signature = 'Assignment',function(assignment){
  if (assignment@log$verbose == T) {
    startTime <- proc.time()
    cat(blue('Calculating correlations '),cli::symbol$continue,'\r',sep = '')
  }
  
  p <- analysisParameters('correlations')
  p@correlations <- assignment@parameters@correlations
  cors <- metabolyse(assignment@data,tibble(ID = 1:nrow(assignment@data)),p,verbose = F) %>%
    correlationResults()
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
    cat(blue('Calculating correlations '),'\t',green(cli::symbol$tick),' ',ncors,' ',elapsed,'\n',sep = '')
  }
  
  return(assignment)
})