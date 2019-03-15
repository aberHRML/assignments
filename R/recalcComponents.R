#' @importFrom tidygraph to_components
#' @importFrom dplyr n

recalcComponents <- function(graph){
  g <- graph %>%
    activate(nodes) %>%
    mutate(Component = group_components()) 
  
  comp <- g %>%
    nodes() %>%
    .$Component %>%
    unique()
  
  weights <- comp %>%
    map(~{
      component <- .
      g %>%
        filter(Component == component) %>%
        edges() %>% 
        .$r %>% 
        mean() %>%
        tibble(Weight = .)
    }) %>%
    set_names(comp) %>%
    bind_rows(.id = 'Component') %>%
    mutate(Component = as.numeric(Component))
  
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