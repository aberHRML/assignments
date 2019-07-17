#' @importFrom metabolyseR analysisParameters metabolyse correlationResults

setMethod('calcCorrelations',signature = 'Assignment',function(assignment){
  p <- analysisParameters('correlations')
  p@correlations <- assignment@parameters@correlations
  cors <- metabolyse(assignment@data,tibble(ID = 1:nrow(assignment@data)),p,verbose = F) %>%
    correlationResults()
  assignment@correlations <- cors
  return(assignment)
})