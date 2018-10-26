#' @importFrom tidygraph as_tbl_graph activate morph unmorph graph_size group_components

calcComponents <- function(MFs,rel) {
  nodes <- MFs
  
  edges <- rel
  
  graph <- as_tbl_graph(edges,directed = F) %>%
    activate(nodes) %>%
    left_join(nodes,by = c('name' = 'Name')) %>%
    mutate(Component = group_components()) %>%
    morph(to_components) %>%
    mutate(Size = graph_size(),
           Nodes = n(),
           Density = (2 * Size) / (Nodes * (Nodes - 1)),
           AverageAddIsoScore = sum(AddIsoScore) / Nodes,
           Plausibility = AverageAddIsoScore * Size) %>%
    unmorph()
  
  
  return(graph)
}