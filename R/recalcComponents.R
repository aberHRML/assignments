#' @importFrom tidygraph to_components
#' @importFrom dplyr n

recalcComponents <- function(graph,parameters){
  g <- graph %>%
    activate(nodes)
  
  comp <- g %>%
    nodes() %>%
    .$Component %>%
    unique()
  
  slaves <- length(comp) / 500
  slaves <-  ceiling(slaves)
  
  if (slaves > parameters@nCores) {
    slaves <- parameters@nCores
  }
  
  clus <- makeCluster(slaves,type = parameters@clusterType)
  
  weights <- comp %>%
    parLapply(cl = clus,function(component,graph){
      graph %>%
        filter(Component == component) %>%
        edges() %>% 
        .$r %>% 
        mean() %>%
        tibble(Weight = .)
    },graph = g) %>%
    set_names(comp) %>%
    bind_rows(.id = 'Component') %>%
    mutate(Component = as.numeric(Component))
  
  stopCluster(clus)
  
  g %>%
    select(-Weight) %>%
    left_join(weights,by = 'Component') %>%
    morph(to_components) %>%
    mutate(Size = graph_size(),
           Nodes = n(),
           Density = (2 * Size) / (Nodes * (Nodes - 1)),
           Weight = sum(Weight) / Nodes,
           AIS = sum(AddIsoScore) / Nodes,
           Plausibility = AIS * Size * Weight) %>%
    unmorph()
}