
plotSolutions <- function(graph,selectedComp,feature){
  graph %>%
    map(~{
      stats <- nodes(.) %>%
        select(Component:Plausibility) %>%
        .[1,]
      
      if (stats$Component[1] == selectedComp){
        border <- 'red'  
      } else {
        border <- 'black'
      }
      
      g <- .
      g <- g %>%
        mutate(Feat = Feature == feature) %>%
        create_layout('nicely')
      ggraph(g) +
        geom_edge_link(aes(colour = r)) +
        scale_edge_color_gradient(low = 'white',high = 'black',limits = c(0.5,1)) +
        geom_node_label(aes(label = name,fill = Feat),size = 2,) +
        scale_fill_manual(values = c('white','steelblue')) +
        theme_graph(base_family = '',
                    title_size = 12,
                    title_face = 'plain',
                    foreground = border,
                    plot_margin = margin(5, 5, 5, 5)) +
        labs(title = str_c('Component ',stats$Component),
             caption = str_c('Degree = ',stats$Degree %>% round(2),'; ',
                             'Weight = ',stats$Weight %>% round(2),'; ',
                             'AIS = ',stats$AIS %>% round(2),'; ',
                             'Plausibility = ',stats$Plausibility %>% round(2))) +
        xlim(min(g$x) - (max(g$x) - min(g$x)) * 0.05,
             max(g$x) + (max(g$x) - min(g$x)) * 0.05) +
        ylim(min(g$y) - (max(g$y) - min(g$y)) * 0.05,
             max(g$y) + (max(g$y) - min(g$y)) * 0.05) +
        guides(fill = 'none')
    }) %>%
    wrap_plots() + plot_annotation(title = str_c('Solutions for feature ',feature))
}

#' plotFeatureSolutions
#' @rdname plotFeatureSolutions
#' @description Plot possible MF solutions for a given feature.
#' @param assignment S4 object of class Assignent
#' @param feature name of feature to plot
#' @param maxComponents maximum number of components to plot
#' @importFrom patchwork plot_annotation
#' @importFrom ggraph create_layout scale_edge_color_gradient geom_node_label
#' @importFrom ggplot2 scale_fill_manual margin xlim ylim guides
#' @export

setMethod('plotFeatureSolutions',signature = 'Assignment',
          function(assignment,feature,maxComponents = 10){
            
            n <- nodes(assignment@addIsoAssign$graph)
            
            comp <- n %>%
              filter(Feature == feature) %>%
              select(Component,Plausibility) %>%
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
              arrange(desc(Plausibility)) %>%
              .$Component
            
            graph <- graph %>%
              set_names(graphComponents) %>%
              .[comp %>% as.character()]
            
            if (length(comp) > maxComponents) {
              graph <- graph[1:maxComponents]
            }
            
            selectedComp <- assignment@addIsoAssign$filteredGraph %>%
              nodes() %>%
              select(Feature,Component) %>%
              filter(Feature == feature) %>%
              .$Component
            
            pl <- plotSolutions(graph,selectedComp,feature)
            
            return(pl)
          })
