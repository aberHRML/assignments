#' @importFrom metabolyseR analysisParameters metabolyse correlationResults keepVariables analysisData dat sinfo
#' @importFrom magrittr set_rownames

setMethod('calcCorrelations',signature = 'Assignment',function(assignment){
  if (assignment@log$verbose == T) {
    startTime <- proc.time()
    cat(blue('Calculating correlations '),cli::symbol$continue,'\r',sep = '')
  }
  
  p <- analysisParameters('correlations')
  p@correlations <- assignment@parameters@correlations
  
  if (str_detect(assignment@parameters@technique,'LC')) {
    feat <- tibble(Feature = colnames(assignment@data)) %>%
      mutate(RT = str_split_fixed(Feature,'@',2)[,2] %>%
               as.numeric())
    RTgroups <- feat %>%
      data.frame() %>%
      set_rownames(.$Feature) %>%
      select(-Feature) %>%
      dist() %>%
      hclust() %>%
      cutree(h = assignment@parameters@RTwindow) %>%
      {tibble(Feature = names(.),Group = .)}
    
    RTsum <- RTgroups %>%
      group_by(Group) %>%
      summarise(N = n())
    
    RTgroups <- RTgroups %>%
      filter(Group %in% {RTsum %>%
          filter(N > 1) %>%
          .$Group})
    
    clus <- makeCluster(assignment@parameters@nCores,type = assignment@parameters@clusterType)
    cors <- RTgroups %>%
      split(.$Group) %>%
      parLapply(cl = clus,function(f){
        analysisData(assignment@data,tibble(ID = 1:nrow(assignment@data))) %>%
          keepVariables(variables = f$Feature) %>%
          {metabolyse(dat(.),sinfo(.),p,verbose = FALSE)} %>% 
          correlationResults()
      }) %>%
      bind_rows(.id = 'RT Group')
    stopCluster(clus)
    
  } else {
    cors <- metabolyse(assignment@data,tibble(ID = 1:nrow(assignment@data)),p,verbose = F) %>%
      correlationResults()  
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
    cat(blue('Calculating correlations '),'\t',green(cli::symbol$tick),' ',ncors,' ',elapsed,'\n',sep = '')
  }
  
  return(assignment)
})