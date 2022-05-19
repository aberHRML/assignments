
graphTheme <- function(){
  ggplot2::theme(
    legend.title = ggplot2::element_text(face = 'bold'),
    plot.margin = ggplot2::margin(5, 5, 5, 5),
    plot.title = ggplot2::element_text(face = 'bold',hjust = 0.5),
    plot.caption = ggtext::element_markdown(hjust = 0.5)
  )
}

plotGraph <- function(graph,min_coef,label_size = 3,axis_offset = 0.1){
  
  g <- graph %>% 
    activate(nodes) %>% 
    mutate(name = str_replace_all(name,'  ','\n') %>% 
             str_replace_all(' ','\n')) %>% 
    ggraph::create_layout('nicely')
  
  g %>% 
    ggraph::ggraph() +
    ggraph::geom_edge_link(ggplot2::aes(colour = coefficient)) +
    ggraph::scale_edge_color_gradient(low = 'lightgrey',
                                      high = 'black',
                                      limits = c(min_coef,1)) +
    ggraph::geom_node_label(ggplot2::aes(label = name),
                            size = label_size) +
    ggraph::theme_graph(base_family = '',
                        base_size = 10) +
    graphTheme() +
    ggplot2::lims(
      x = c(
        min(g$x) - (max(g$x) - min(g$x)) * axis_offset,
        max(g$x) + (max(g$x) - min(g$x)) * axis_offset
      ),
      y = c(
        min(g$y) - (max(g$y) - min(g$y)) * axis_offset,
        max(g$y) + (max(g$y) - min(g$y)) * axis_offset
      )
    )
}

#' Plot a component
#' @rdname plotComponent
#' @description Plot a molecular formula component graph.
#' @param assignment S4 object of class Assignment
#' @param component component number to extract
#' @param iteration the assignment iteration
#' @param type the graph type to return. `filtered` returns the assignment graph after component selection. `all` returns all assignment components.
#' @param label_size node label size
#' @param axis_offset axis proportion by which to increase axis limits. Prevents cut off of node labels.
#' @examples 
#' \dontrun{
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' assignment <- assignMFs(feature_data,p)
#' 
#' plotComponent(assignment,1,'A&I1')
#' }
#' @export

setGeneric('plotComponent',
           function(assignment,
                    component,
                    iteration,
                    type = c('filtered','all'),
                    label_size = 3,
                    axis_offset = 0.1)
             standardGeneric('plotComponent'))

#' @importFrom dplyr mutate_if

setMethod('plotComponent',signature = 'Assignment',
          function(assignment,
                   component,
                   iteration,
                   type = c('filtered','all'),
                   label_size = 3,
                   axis_offset = 0.1
          ){
            
            check_installed(c('ggraph',
                              'ggplot2',
                              'ggtext'))
            
            component_graph <- component(assignment,
                                         component,
                                         iteration,
                                         type) 
            
            component_stats <- component_graph %>% 
              nodes() %>% 
              select(`MF Plausibility (%)`,
                     AIS,
                     Component:`Component Plausibility`) %>% 
              distinct() %>% 
              mutate_if(is.numeric,signif,digits = 3)
            
            min_coef <- correlationsParameters(assignment)$minCoef
            
            plotGraph(component_graph,
                      min_coef,
                      label_size,
                      axis_offset) +
              ggplot2::labs(
                title = paste0('Component ',component),
                caption =  glue::glue('
          P<sub>c</sub> = {component_stats$`Component Plausibility`};
          Degree = {component_stats$Degree};
          
          AIS<sub>c</sub> = {component_stats$AIS};
          P<sub>MF</sub> = {component_stats$`MF Plausibility (%)`}%')
            )
          })

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
