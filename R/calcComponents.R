#' @importFrom tidygraph as_tbl_graph activate morph unmorph graph_size group_components
#' @importFrom magrittr set_names

calcComponents <- function(MFs,rel,nCores,clusterType) {
  no <- MFs
  
  ed <- rel
  
  graph <- as_tbl_graph(ed,directed = F) %>%
    activate(nodes) %>%
    left_join(no,by = c('name' = 'Name')) %>%
    mutate(Component = group_components()) 
  
  comp <- graph %>%
    nodes() %>%
    .$Component %>%
    unique()
  
  slaves <- {length(comp) / 1000} %>%
    ceiling()
  
  if (slaves > nCores) {
    slaves <- nCores
  }
  
  clus <- makeCluster(slaves, type = clusterType)
  
  weights <- comp %>%
    parLapply(cl = clus,fun = function(x,graph){
      component <- x
      graph %>%
        filter(Component == component) %>%
        edges() %>% 
        .$r %>% 
        mean() %>%
        tibble(Weight = .)
      },graph = graph) %>%
    set_names(comp) %>%
    bind_rows(.id = 'Component') %>%
    mutate(Component = as.numeric(Component))
  
  stopCluster(clus)
  
  graph <- graph %>%
    left_join(weights,by = 'Component') %>%
    morph(to_components) %>%
    mutate(Size = graph_size(),
           Nodes = n(),
           Density = (2 * Size) / (Nodes * (Nodes - 1)),
           Weight = sum(Weight) / Nodes,
           AIS = sum(AddIsoScore) / Nodes,
           Plausibility = AIS * Size * Weight) %>%
    unmorph()
  
  return(graph)
}