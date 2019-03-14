#' @importFrom tidygraph to_components
#' @importFrom dplyr n

recalcComponents <- function(graph){
  graph %>%
    activate(nodes) %>%
    morph(to_components) %>%
    mutate(Size = graph_size(),
           Nodes = n(),
           Density = (2 * Size) / (Nodes * (Nodes - 1)),
           AverageAddIsoScore = sum(AddIsoScore) / Nodes,
           Plausibility = AverageAddIsoScore * Size) %>%
    unmorph() %>%
    filter(Nodes > 1)
}