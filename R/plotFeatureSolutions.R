#' plotFeatureSolutions
#' @rdname plotFeatureSolutions
#' @importFrom patchwork plot_annotation
#' @export

setMethod('plotFeatureSolutions',signature = 'Assignment',
          function(assignment,feature){
            
            n <- nodes(assignment@addIsoAssign$graph)
            
            comp <- n %>%
              filter(Feature == feature) %>%
              .$Component %>%
              unique()
            
            graph <- assignment@addIsoAssign$graph %>%
              filter(Component %in% comp) %>%
              mutate(name = str_replace_all(name,'  ','\n')) %>%
              mutate(name = str_replace_all(name,' ','\n')) %>%
              morph(to_components)
            
            selectedComp <- assignment@addIsoAssign$filteredGraph %>%
              nodes() %>%
              .$Component %>%
              unique() %>%
              {comp[comp %in% .]}
            
            graph %>%
              map(~{
                stats <- nodes(.) %>%
                  select(Component:Plausibility) %>%
                  .[1,]
                
                if (stats$Component == selectedComp){
                  border <- 'red'  
                } else {
                  border <- 'black'
                }
                
                g <- .
                g <- g %>%
                  mutate(Feat = Feature == feature) %>%
                  create_layout('nicely')
                ggraph(g) +
                  geom_edge_link() +
                  geom_node_label(aes(label = name,colour = Feat),size = 2,) +
                  theme_graph(title_size = 12,
                              title_face = 'plain',
                              foreground = border,
                              plot_margin = margin(5, 5, 5, 5)) +
                  labs(title = str_c('Component ',stats$Component),
                       caption = str_c('Size = ',stats$Size,'; ',
                                       'Nodes = ',stats$Nodes,'; ',
                                       'Density = ',stats$Density %>% round(2),'; ',
                                       'AddIsoScore = ',stats$AverageAddIsoScore %>% round(2),'; ',
                                       'Plausibility = ',stats$Plausibility %>% round(2))) +
                  xlim(min(g$x) - (max(g$x) - min(g$x)) * 0.05,
                       max(g$x) + (max(g$x) - min(g$x)) * 0.05) +
                  ylim(min(g$y) - (max(g$y) - min(g$y)) * 0.05,
                       max(g$y) + (max(g$y) - min(g$y)) * 0.05) +
                  guides(colour = FALSE) +
                  scale_colour_manual(values = c('black','steelblue'))
              }) %>%
              .[order(comp)] %>%
              wrap_plots() + plot_annotation(title = str_c('Solutions for feature ',feature))
          })