
networkPlot <- function(network,layout,rThreshold,assignedNodes,explainedEdges){
  check_installed(c('ggraph',
                    'ggthemes',
                    'ggplot2',
                    'graphlayouts'))
  
  rt <- str_c('Visualised using threshold of r > ',rThreshold)
  nn <- str_c('Total nodes = ',sum(assignedNodes))
  an <- str_c('Assigned nodes = ',
              assignedNodes[1],
              ' (',
              {assignedNodes[1]/sum(assignedNodes) * 100} %>%
                round(),'%)')
  ne <- str_c('Total edges = ',sum(explainedEdges))
  ee <- str_c('Explained edges = ',
              explainedEdges[1],
              ' (',
              {explainedEdges[1]/sum(explainedEdges) * 100} %>%
                round(),'%)')
  
  ggraph::ggraph(network,layout = layout) +
    ggraph::geom_edge_link(alpha = 0.2) +
    ggraph::geom_node_point(ggplot2::aes(fill = Assigned),shape = 21) +
    ggthemes::scale_fill_ptol() +
    ggraph::theme_graph() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::coord_fixed() +
    ggraph::facet_edges(~Explained) +
    ggplot2::labs(title = str_c('Assignment correlation network'),
         caption = str_c(rt,nn,an,ne,ee,sep = '\n'))
}

#' Plot an assignment network
#' @rdname plotNetwork
#' @description plot assignment network
#' @param assignment of class Assignment
#' @param layout graph layout to use. See \code{\link[ggraph]{ggraph}} for layout options
#' @param rThreshold r threhold to use for filtering edge correlation weights
#' @importFrom tidygraph as_tbl_graph bind_graphs
#' @importFrom igraph set_vertex_attr set_edge_attr
#' @export

setGeneric('plotNetwork',function(assignment, layout = 'stress', rThreshold = 0.7){
  standardGeneric('plotNetwork')
})

#' @rdname plotNetwork

setMethod('plotNetwork',signature = 'Assignment',
          function(assignment, layout = 'stress', rThreshold = 0.7){
            
            AI <- assignment@addIsoAssign$filtered_graph
            TA <- assignment@transAssign %>%
              map(~{.$filtered_graph})
            
            if (length(TA) > 0){
              if (length(TA) > 1) {
                graph <- bind_graphs(graph,
                  {a <- TA[[1]]
                  for (i in 2:length(TA)) {
                    a <- bind_graphs(a,TA[[i]])
                  }
                  a
                  }
                ) 
              } else {
                graph <- bind_graphs(AI,TA[[1]])
              }
            } 
            
            e <- edges(graph) %>%
              mutate(Explained = 'Explained')
            n <- nodes(graph) %>%
              select(name:`MF Plausibility (%)`) %>%
              distinct() %>%
              mutate(Assigned = 'Assigned')
            
            network <- assignment %>%
              .@correlations %>%
              filter(coefficient > rThreshold) %>%
              as_tbl_graph(directed = F) %>%
              activate(nodes) %>%
              rename(Feature = name) %>%
              mutate(Mode = str_sub(Feature,1,1)) %>%
              left_join(n, by = "Feature") %>%
              activate(edges) %>%
              left_join(e, by = c("Mode1", "Mode2", "m/z1", "m/z2", "RetentionTime1", "RetentionTime2", "log2IntensityRatio", "coefficient", "ID"))
            
            assigned <- nodes(network)$Assigned
            assigned[is.na(assigned)] <- 'Unassigned'
            
            network <- set_vertex_attr(network,'Assigned',value = assigned)
            
            explained <- edges(network)$Explained
            explained[is.na(explained)] <- 'Unexplained'
            
            network <- set_edge_attr(network,'Explained',value = explained)
            
            explainedEdges <- network %>% 
              edges() %>%
              .$Explained %>%
              table()
            
            assignedNodes <- network %>% 
              nodes() %>%
              .$Assigned %>%
              table()
            
            networkPlot(network,
                        layout,
                        rThreshold,
                        assignedNodes,
                        explainedEdges)
          })
            