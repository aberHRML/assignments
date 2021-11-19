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