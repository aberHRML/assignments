
degree <- function(n_nodes,n_edges){
  2 * (n_edges / n_nodes)
}

plausibility <- function(AIS,degree,weight){
  AIS + weight
}

componentMetrics <- function(component,max_add_iso_total){
  component %>% 
  mutate(Size = graph_size(),
         Nodes = n(),
         Degree = degree(Nodes,Size),,
         Density = (2 * Size) / (Nodes * (Nodes - 1)),
         Weight = sum(Weight) / Nodes,
         AIS = Size * sum(AddIsoScore) / max_add_iso_total,
         `Component Plausibility` = plausibility(AIS,Degree,Weight))
}

componentFilters <- function(){
  tibble(Measure = c('Component Plausibility',
                     'Degree',
                     'AIS',
                     'MF Plausibility (%)',
                     'PPM error'),
         Direction = c(rep('max',4),'min'))
}

#' @importFrom tidygraph as_tbl_graph activate morph unmorph graph_size group_components tbl_graph
#' @importFrom magrittr set_names

calcComponents <- function(graph_nodes,
                           graph_edges,
                           assignment) {
  
  graph <- as_tbl_graph(graph_edges,directed = FALSE) %>%
    activate(nodes) %>%
    left_join(graph_nodes,by = c('name' = 'Name')) %>%
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
    componentMetrics(max_add_iso_total = maxAddIsoScore(assignment)) %>%
    unmorph()
  
  return(graph)
}

#' @importFrom tidygraph to_components
#' @importFrom dplyr n

recalcComponents <- function(graph,
                             assignment){
  g <- graph %>%
    activate(nodes)
  
  comp <- g %>%
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
    },graph = g) %>%
    set_names(comp) %>%
    bind_rows(.id = 'Component') %>%
    mutate(Component = as.numeric(Component))
  
  g %>%
    select(-Weight) %>%
    left_join(weights,by = 'Component') %>%
    morph(to_components) %>%
    componentMetrics(max_add_iso_total = maxAddIsoScore(assignment)) %>% 
    unmorph()
}

filterComponents <- function(graph,
                             assignment,
                             filters = componentFilters()){
  filtered_graph <- graph
  
  for (i in 1:nrow(filters)) { 
    f <- filters[i,]
    filtered_graph <- filtered_graph %>%
      activate(nodes) %>%
      filter(name %in% {filtered_graph %>% 
          nodes() %>% 
          eliminate(f$Measure,f$Direction) %>%
          .$name}) 
    if (V(filtered_graph) %>% length() > 0) {
      filteredGraph <- filtered_graph %>%
        recalcComponents(assignment)
    } else {
      break()
    }
  }
  
  return(filtered_graph)
}