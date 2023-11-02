
avg_degree <- function(n_nodes,n_edges){
  2 * (n_edges / n_nodes)
}

plausibility <- function(degree,AIS,weight){
  degree * AIS * weight
}

#' @importFrom purrr compact
#' @importFrom tidygraph bind_graphs

clean <- function(graph,adduct_rules_table){
  cleaned_graph <- graph %>% 
    morph(to_components) %>% 
    furrr::future_map(~{
      component_adducts <- .x %>% 
        nodes() %>% 
        .$Adduct %>%
        unique()
      
      adduct_info <- adduct_rules_table %>%
        select(Adduct = Name,
               adduct_isotopic = Isotopic)
      
      component_nodes <- .x %>% 
        nodes() %>% 
        left_join(adduct_info, 
                  by = "Adduct") %>% 
        mutate(isotopic = !is.na(Isotope),
               adduct_isotopic = as.logical(adduct_isotopic),
               either_isotopic = adduct_isotopic | isotopic) 
      
      if (all(component_nodes$isotopic) == TRUE |
          all(component_nodes$adduct_isotopic) == TRUE |
          all(component_nodes$either_isotopic) == TRUE |
          nrow(component_nodes) < 2) NULL
      else return(.x)
    }) %>% 
    compact()
  
  if (length(cleaned_graph) == 0){
    cleaned_graph <- graph %>% 
      slice(0)
  }
  
  if (length(cleaned_graph) > 0){
    cleaned_graph <- cleaned_graph %>% 
      bind_graphs() 
  }
   
  return(cleaned_graph)
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
         Degree = avg_degree(Nodes,Size),
         Density = (2 * Size) / (Nodes * (Nodes - 1)),
         Weight = sum(Weight) / Nodes,
         AIS = sum(AIS) / max_add_iso_total,
         `Component Plausibility` = plausibility(Degree,AIS,Weight)
         )
}

componentFilters <- function(){
  tibble(Measure = c('Component Plausibility',
                     'MF Plausibility (%)',
                     'PPM error'),
         Direction = c(rep('max',2),'min'))
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
    clean(adductRules(assignment))
  
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
          abs() %>% 
          mean() %>%
          tibble(Weight = .)
      },graph = graph) %>%
      set_names(comp) %>%
      bind_rows(.id = 'Component') %>%
      mutate(Component = as.numeric(Component))
    
    graph <- graph %>%
      left_join(weights,by = 'Component') %>%
      morph(to_components) %>%
      componentMetrics(max_add_iso_total = maxAIS(assignment)) %>%
      unmorph() 
  }
  
  return(graph)
}

#' @importFrom tidygraph to_components
#' @importFrom dplyr n

recalcComponents <- function(graph,
                             assignment){
  
  graph <- graph %>% 
    clean(adductRules(assignment))
  
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
      componentMetrics(max_add_iso_total = maxAIS(assignment)) %>% 
      unmorph() 
  } 
  
  return(graph)
}

#' @importFrom igraph degree

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
          .$name}) %>% 
      filter(degree(.) != 0)
    if (E(filtered_graph) %>% length() > 0) {
      filtered_graph <- filtered_graph %>%
        recalcComponents(assignment)
    } else {
      break()
    }
  }
  
  return(filtered_graph)
}
