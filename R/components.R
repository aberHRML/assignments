
degree <- function(n_nodes,n_edges){
  2 * (n_edges / n_nodes)
}

plausibility <- function(AIS,degree,weight){
  AIS + weight
}

#' @importFrom purrr compact

clean <- function(graph){
  graph %>% 
    morph(to_components) %>% 
    map(~{
      component_nodes <- nodes(.x)
      
      if (nrow(component_nodes) > 1 &
          NA %in% component_nodes$Isotope) return(.x)
      else NULL
    }) %>% 
    compact() %>% 
    bind_graphs() 
}

nComponents <- function(graph){
  graph %>% 
    morph(to_components) %>% 
    length()
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
    mutate(Component = group_components()) %>% 
    clean()
  
  if (nComponents(graph) > 0){
    comp <- graph %>%
      nodes() %>%
      .$Component %>%
      unique()
    
    weights <- comp %>%
      future_map(~{
        graph %>%
          filter(Component == .x) %>%
          edges() %>% 
          .$coefficient %>% 
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
  }
  
  return(graph)
}

#' @importFrom tidygraph to_components
#' @importFrom dplyr n

recalcComponents <- function(graph,
                             assignment){
  
  graph <- graph %>% 
    clean()
  
  if (nComponents(graph) > 0){
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
          .$coefficient %>% 
          mean() %>%
          tibble(Weight = .)
      },graph = g) %>%
      set_names(comp) %>%
      bind_rows(.id = 'Component') %>%
      mutate(Component = as.numeric(Component))
    
    graph <- g %>%
      select(-Weight) %>%
      left_join(weights,by = 'Component') %>%
      morph(to_components) %>%
      componentMetrics(max_add_iso_total = maxAddIsoScore(assignment)) %>% 
      unmorph() 
  } 
  
  return(graph)
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
      filtered_graph <- filtered_graph %>%
        recalcComponents(assignment)
    } else {
      break()
    }
  }
  
  return(filtered_graph)
}