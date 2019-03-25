#' @importFrom tidygraph to_components
#' @importFrom dplyr n

recalcComponents <- function(graph,nCores,clusterType){
  g <- graph %>%
    activate(nodes)
  
  comp <- g %>%
    nodes() %>%
    .$Component %>%
    unique()
  
  slaves <- {length(comp) / 1000} %>%
    ceiling()
  
  if (slaves > nCores) {
    slaves <- nCores
  }
  
  clus <- makeCluster(slaves, type = clusterType)
  clusterExport(clus,c('filter','edges','tibble'))
  
  weights <- comp %>%
    parLapply(cl = clus,fun = function(x,graph){
      component <- x
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