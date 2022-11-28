
graphTheme <- function(){
  ggplot2::theme(
    legend.title = ggplot2::element_text(face = 'bold'),
    plot.margin = ggplot2::margin(5, 5, 5, 5),
    plot.title = ggplot2::element_text(face = 'bold',hjust = 0.5),
    plot.caption = ggtext::element_markdown(hjust = 0.5)
  )
}

plotGraph <- function(graph,
                      min_coef,
                      label_size = 3,
                      axis_offset = 0.1,
                      border = 'black',
                      highlight = NA){
  
  if (!is.na(highlight)){
    graph <- graph %>% 
      activate(nodes) %>% 
      mutate(selected = Feature == highlight)
  }
  
  g <- graph %>% 
    activate(nodes) %>% 
    mutate(name = str_replace_all(name,'  ','\n') %>% 
             str_replace_all(' ','\n')) %>% 
    ggraph::create_layout('nicely')
  
  pl <- g %>% 
    ggraph::ggraph() +
    ggraph::geom_edge_link(ggplot2::aes(colour = coefficient)) +
    ggraph::scale_edge_color_gradient(low = 'lightgrey',
                                      high = 'black',
                                      limits = c(min_coef,1))
  
  if (!is.na(highlight)) {
    pl <- pl + 
      ggraph::geom_node_label(
        ggplot2::aes(label = name,fill = selected),
        size = label_size) +
      ggplot2::scale_fill_manual(values = c('white','lightblue')) +
      ggplot2::guides(fill = 'none')
  } else {
    pl <- pl + 
      ggraph::geom_node_label(
        ggplot2::aes(label = name),
        size = label_size)
  }
  
  pl + 
    ggraph::theme_graph(base_family = '',
                        base_size = 10,
                        title_size = 11,
                        foreground = border) +
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
#' @param type the graph type to return. `selected` returns the assignment graph after component selection. `all` returns all assignment components.
#' @param label_size node label size
#' @param axis_offset axis proportion by which to increase axis limits. Prevents cut off of node labels.
#' @param border specify a plot border colour
#' @param highlight specify a feature node to highlight
#' @examples 
#' \dontrun{
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' mf_assignments <- assignMFs(feature_data,p)
#' 
#' plotComponent(mf_assignments,1,'A&I1')
#' }
#' @export

setGeneric('plotComponent',
           function(assignment,
                    component,
                    iteration,
                    type = c('selected','all'),
                    label_size = 3,
                    axis_offset = 0.1,
                    border = NA,
                    highlight = NA)
             standardGeneric('plotComponent'))

#' @importFrom dplyr mutate_if
#' @rdname plotComponent

setMethod('plotComponent',signature = 'Assignment',
          function(assignment,
                   component,
                   iteration,
                   type = c('selected','all'),
                   label_size = 3,
                   axis_offset = 0.1,
                   border = NA,
                   highlight = NA 
          ){
            
            check_installed(c('ggraph',
                              'ggplot2',
                              'ggtext',
                              'glue'))
            
            component_graph <- component(assignment,
                                         component,
                                         iteration,
                                         type) 
            
            if (!is.na(highlight) & 
                !highlight %in% nodes(component_graph)$Feature) {
              stop(paste0('Highlight feature ',highlight,' not found in component.'))
            }
            
            component_stats <- component_graph %>% 
              nodes() %>% 
              select(AIS,
                     Component:`Component Plausibility`) %>% 
              distinct() %>% 
              mutate_if(is.numeric,signif,digits = 3)
            
            min_coef <- correlationsParameters(assignment)$minCoef
            
            plotGraph(component_graph,
                      min_coef,
                      label_size,
                      axis_offset,
                      border,
                      highlight
                      ) +
              ggplot2::labs(
                title = paste0('Component ',component),
                caption =  glue::glue('
          P<sub>c</sub> = {component_stats$`Component Plausibility`};
          Degree = {component_stats$Degree};
          AIS<sub>c</sub> = {component_stats$AIS}'
              ))
          })

#' Plot the solutions for a feature
#' @rdname plotFeatureComponents
#' @description Plot possible MF solutions for a given feature.
#' @param assignment S4 object of class Assignent
#' @param feature name of feature to plot
#' @param iteration components from which iteration to plot
#' @param type the graph type to return. `all` returns all assignment components. `selected` returns the assignment graph after component selection.
#' @param max_components maximum number of components to plot
#' @param label_size node label size
#' @param axis_offset axis proportion by which to increase axis limits. Prevents cut off of node labels.
#' @export

setGeneric('plotFeatureComponents',
           function(assignment,
                    feature,
                    iteration,
                    type = c('all','selected'),
                    max_components = 6,
                    label_size = 3,
                    axis_offset = 0.1)
             standardGeneric('plotFeatureComponents')
)

#' @rdname plotFeatureComponents
#' @importFrom dplyr slice

setMethod('plotFeatureComponents',signature = 'Assignment',
          function(assignment,
                   feature,
                   iteration,
                   type = c('all','selected'),
                   max_components = 6,
                   label_size = 2,
                   axis_offset = 0.05){
            
            check_installed(c('ggraph',
                              'ggplot2',
                              'ggtext',
                              'patchwork'))
            
            if (!feature %in% colnames(featureData(assignment))) {
              stop('Feature not found in assignment data.',
                   call. = FALSE)
            }
            
            type <- match.arg(type,
                              choices = c('all','selected'))
            
            selected_component <- assignments(assignment) %>% 
              filter(Feature == feature,
                     Iteration == iteration) %>% 
              .$Component
            
            feature_components <- featureComponents(assignment,feature,type) %>% 
              filter(Iteration == iteration) %>% 
              select(Component) %>% 
              arrange(Component) %>% 
              mutate(border = 'black',
                     border = border %>% 
                       replace(Component == selected_component,
                               'red'))
            
            if (nrow(feature_components) == 0){
              stop(paste0('No components for feature ',
                          feature,
                          ' found in iteration ',
                          iteration,'.'),
                   call. = FALSE)
            }
            
            if (nrow(feature_components) > max_components){
              feature_components <- slice(
                feature_components,
                seq_len(max_components))
            }
            
            pl <- feature_components %>%
              rowwise() %>% 
              group_split() %>% 
              map(~plotComponent(
                assignment,
                .x$Component,
                iteration,
                type,
                label_size,
                axis_offset,
                highlight = feature,
                border = .x$border
              )) %>% 
              patchwork::wrap_plots() +
              patchwork::plot_layout(guides = 'collect')
            
            if (length(selected_component) > 0){
              pl <- pl +
                patchwork::plot_annotation(
                  caption = 'Red highlighted graph denotes the component selected for assignment.',
                  theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
                )
            }
            
            return(pl)
          })
