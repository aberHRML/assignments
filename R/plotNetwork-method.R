#' plotNetwork-Assignment
#' @rdname plotNetwork
#' @description plot assignment network
#' @param assignment of class Assignment
#' @param layout graph layout to use. See \code{\link[ggraph]{ggraph}} for layout options
#' @param rThreshold r threhold to use for filtering edge correlation weights
#' @importFrom tidygraph as_tbl_graph bind_graphs
#' @importFrom igraph set_vertex_attr set_edge_attr
#' @importFrom ggraph ggraph geom_edge_link geom_node_point theme_graph geom_node_text facet_edges
#' @importFrom ggthemes scale_fill_ptol
#' @importFrom ggplot2 labs aes element_blank coord_fixed
#' @importFrom graphlayouts layout_igraph_stress
#' @export

setMethod('plotNetwork',signature = 'Assignment',
          function(assignment, layout = 'stress', rThreshold = 0.7){
            
            AI <- assignment@addIsoAssign$filteredGraph
            TA <- assignment@transAssign %>%
              map(~{.$filteredGraph})
            
            graph <- AI %>%
              bind_graphs({a <- TA[[1]]
              for (i in 2:length(TA)) {
                a <- bind_graphs(a,TA[[i]])
              }
              a
              }) 
            
            e <- edges(graph) %>%
              mutate(Explained = 'Explained')
            n <- nodes(graph) %>%
              select(name:Score) %>%
              distinct() %>%
              mutate(Assigned = 'Assigned')
            
            network <- assignment %>%
              .@preparedCorrelations %>%
              filter(r > rThreshold) %>%
              as_tbl_graph(directed = F) %>%
              activate(nodes) %>%
              rename(Feature = name) %>%
              mutate(Mode = str_sub(Feature,1,1)) %>%
              left_join(n, by = "Feature") %>%
              activate(edges) %>%
              left_join(e, by = c("Mode1", "Mode2", "m/z1", "m/z2", "RetentionTime1", "RetentionTime2", "log2IntensityRatio", "r", "ID"))
            
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
            
            ggraph(network,layout = layout) +
              geom_edge_link(alpha = 0.2) +
              geom_node_point(aes(fill = Assigned),shape = 21) +
              scale_fill_ptol() +
              theme_graph() +
              theme(legend.title = element_blank()) +
              coord_fixed() +
              facet_edges(~Explained) +
              labs(title = str_c('Assignment correlation network'),
                   caption = str_c(rt,nn,an,ne,ee,sep = '\n'))
          })