#' @importFrom tidygraph as_tbl_graph activate morph unmorph graph_size group_components tbl_graph
#' @importFrom magrittr set_names

calcComponents <- function(MFs,rel,parameters) {
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
  
  weights <- comp %>%
    future_map(~{
      graph %>%
        filter(Component == .x) %>%
        edges() %>% 
        .$r %>% 
        mean() %>%
        tibble(Weight = .)
      },graph = graph) %>%
    set_names(comp) %>%
    bind_rows(.id = 'Component') %>%
    mutate(Component = as.numeric(Component))
  
  graph <- graph %>%
    left_join(weights,by = 'Component') %>%
    morph(to_components) %>%
    mutate(Size = graph_size(),
           Nodes = n(),
           Density = (2 * Size) / (Nodes * (Nodes - 1)),
           Weight = sum(Weight) / Nodes,
           AIS = sum(AddIsoScore) / Nodes,
           Plausibility = plausibility(AIS,Size,Weight)) %>%
    unmorph()
  
  return(graph)
}