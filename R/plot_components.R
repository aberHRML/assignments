
plotSolutions <- function(graph,selectedComp,feature){
  check_installed(c('ggraph',
                    'ggplot2',
                    'patchwork'))
  graph %>%
    map(~{
      stats <- nodes(.x) %>%
        select(Component:`Component Plausibility`,AIS) %>%
        .[1,]
      
      if (length(selectedComp) > 0){
        border <- ifelse(stats$Component[1] == selectedComp,
                         'red',
                         'black')
      } else {
        border <- 'black'
      }
      
      g <- .x %>%
        mutate(Feat = Feature == feature) %>%
        ggraph::create_layout('nicely')
      ggraph::ggraph(g) +
        ggraph::geom_edge_link(ggplot2::aes(colour = coefficient)) +
        ggraph::scale_edge_color_gradient(low = 'white',
                                          high = 'black',
                                          limits = c(0.5,1)) +
        ggraph::geom_node_label(ggplot2::aes(label = name,
                                             fill = Feat),
                                size = 2,) +
        ggplot2::scale_fill_manual(values = c('white','steelblue')) +
        ggraph::theme_graph(base_family = '',
                            title_size = 12,
                            title_face = 'plain',
                            foreground = border,
                            plot_margin = ggplot2::margin(5, 5, 5, 5)) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = 'bold',
                                             hjust = 0.5),
          plot.caption = ggplot2::element_text(hjust = 0)) +
        ggplot2::labs(
          title = str_c('Component ',stats$Component),
          caption = str_c('Degree = ',stats$Degree %>% round(2),'; ',
                          'Weight = ',stats$Weight %>% round(2),'; ',
                          'AIS = ',stats$AIS %>% round(2),'; ',
                          'Plausibility = ',stats$`Component Plausibility` %>% 
                            round(2))) +
        ggplot2::xlim(min(g$x) - (max(g$x) - min(g$x)) * 0.05,
                      max(g$x) + (max(g$x) - min(g$x)) * 0.05) +
        ggplot2::ylim(min(g$y) - (max(g$y) - min(g$y)) * 0.05,
                      max(g$y) + (max(g$y) - min(g$y)) * 0.05) +
        ggplot2::guides(fill = 'none')
    }) %>%
    patchwork::wrap_plots() + 
    patchwork::plot_layout(guides = 'collect') +
    patchwork::plot_annotation(
      title = str_c('Solutions for feature ',
                    feature),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = 'bold',
                                           hjust = 0.5)))
}

#' Plot the solutions for a feature
#' @rdname plotFeatureSolutions
#' @description Plot possible MF solutions for a given feature.
#' @param assignment S4 object of class Assignent
#' @param feature name of feature to plot
#' @param maxComponents maximum number of components to plot
#' @export

setGeneric('plotFeatureSolutions',
           function(assignment,feature,maxComponents = 10)
             standardGeneric('plotFeatureSolutions')
)

#' @rdname plotFeatureSolutions

setMethod('plotFeatureSolutions',signature = 'Assignment',
          function(assignment,feature,maxComponents = 10){
            
            if (!feature %in% colnames(featureData(assignment))) {
              stop('Feature not found in assignment data')
            }
            
            available_iterations <- iterations(assignment)
            
            graphs <- available_iterations %>% 
              map(graph,assignment = assignment,type = 'all') %>% 
              set_names(available_iterations)
            
            available_nodes <- graphs %>% 
              map_dfr(nodes,.id = 'Iteration') %>% 
              select(Feature) %>% 
              distinct()
            
            if (!(feature %in% available_nodes$Feature)){
              stop(
                paste0('No assignment solutions available for feature ',
                       feature),
                   call. = FALSE)
            }
            
            comp <- n %>%
              filter(Feature == feature) %>%
              select(Component,`Component Plausibility`) %>%
              distinct() %>%
              arrange(Component)
            
            graph <- assignment@addIsoAssign$graph %>%
              filter(Component %in% comp$Component) %>%
              mutate(name = str_replace_all(name,'  ','\n')) %>%
              mutate(name = str_replace_all(name,' ','\n')) %>%
              morph(to_components) 
            
            graphComponents <- graph %>% 
              map_dbl(~{nodes(.) %>% 
                  .$Component %>% 
                  .[1]
              })
            
            comp <- comp %>%
              arrange(desc(`Component Plausibility`)) %>%
              .$Component
            
            graph <- graph %>%
              set_names(graphComponents) %>%
              .[comp %>% as.character()]
            
            if (length(comp) > maxComponents) {
              graph <- graph[1:maxComponents]
            }
            
            selectedComp <- assignment@addIsoAssign$filtered_graph %>%
              nodes() %>%
              select(Feature,Component) %>%
              filter(Feature == feature) %>%
              .$Component
            
            pl <- plotSolutions(graph,selectedComp,feature)
            
            return(pl)
          })
