#' @importFrom stringr str_c

setGeneric("transformationAssign", function(assignment)
  standardGeneric("transformationAssign"))

#' @importFrom dplyr full_join select distinct group_split
#' @importFrom mzAnnotation transformMF

setMethod('transformationAssign',signature = 'Assignment',
          function(assignment){
            
            if (assignment@log$verbose == TRUE) {
              t_start_time <- proc.time()
              message(blue('Transformation assignment'),
                      cli::symbol$continue)
            }
            count <- 0
            repeat {
              count <- count + 1
              
              if (assignment@log$verbose == TRUE) {
                start_time <- proc.time()
                message(str_c('iteration ', 
                              count,' '),
                        cli::symbol$continue,
                        '\r',
                        appendLF = FALSE)
              }
              
              assigned <- assignments(assignment)
              
              rel <- assignment %>% 
                relationships() %>%
                filter(
                  (`m/z1` %in% assigned$`Measured m/z` | 
                     (`m/z2` %in% assigned$`Measured m/z`)) & 
                    !(`m/z1` %in% assigned$`Measured m/z` & 
                        (`m/z2` %in% assigned$`Measured m/z`)),
                  !(is.na(Transformation1) & is.na(Transformation2))
                )
              
              mz1 <- rel %>%
                semi_join(assigned,
                          by = c('m/z1' = 'Measured m/z', 
                                 'Adduct1' = 'Adduct', 
                                 'Isotope1' = 'Isotope')) %>%
                filter(is.na(Transformation1))
              mz2 <- rel %>%
                semi_join(assigned,
                          by = c('m/z2' = 'Measured m/z',
                                 'Adduct2' = 'Adduct',
                                 'Isotope2' = 'Isotope')) %>%
                filter(is.na(Transformation2))
              
              rel <- bind_rows(mz1,mz2)
              
              if (nrow(rel) == 0) break
              
              M <- collateM(rel,
                            maxM(assignment))%>%
                filter(!(mz %in% assigned$`Measured m/z`))
              
              MFs <- generateMFs(M,
                                 ppm(assignment),
                                 MFrankThreshold(assignment),
                                 adductRules(assignment),
                                 isotopeRules(assignment),
                                 AIS(assignment))
              
              if (nrow(MFs) == 0) break
              
              MFs <- MFs %>%
                bind_rows(assigned %>% 
                            select(dplyr::any_of(names(MFs))) %>% 
                            left_join(AIS(assignment),
                                      by = c('Adduct','Isotope')))
              graph_edges <- rel %>% 
                addMFs(MFs,
                       identMF = FALSE) %>%
                sanitiseTransformations(assignment@transformation_rules) %>% 
                mutate(RetentionTime1 = as.numeric(RetentionTime1),
                       RetentionTime2 = as.numeric(RetentionTime2)) %>%
                addNames()
              
              if (nrow(graph_edges) == 0) break
              
              graph_nodes <- collateMFs(graph_edges,MFs)
              
              graph <- calcComponents(graph_nodes,
                                      graph_edges,
                                      assignment)
              
              filtered_graph <- filterComponents(graph,
                                                 assignment)
              
              newly_assigned <- filtered_graph %>%
                nodes() %>% 
                rename(Name = name) %>%
                mutate(Mode = str_sub(Feature,1,1)) %>%
                filter(!(Name %in% assigned$Name)) %>%
                select(Name:`MF Plausibility (%)`,Mode) %>%
                mutate(Iteration = str_c('T',count + 1)) %>% 
                group_split(MF) %>% 
                map_dfr(~{
                  if (NA %in% .x$Isotope) return(.x)
                  else NULL
                })
              
              if (nrow(newly_assigned) == 0) break
              
              assignment@transAssign[[count]] <- list(
                graph = graph,
                filtered_graph = filtered_graph,
                assigned = newly_assigned)
              
              assignment@assignments <- bind_rows(assignment@assignments,
                                                  newly_assigned)
              
              if (isTRUE(assignment@log$verbose)) {
                end_time <- proc.time()
                elapsed <- elapsedTime(start_time,end_time)
                message(str_c('iteration ', 
                              count,' '),
                        '\t\t\t',
                        green(cli::symbol$tick),
                        ' ',
                        elapsed)
              }
            }
            
            names(assignment@transAssign) <- paste0('T',
                                                    seq_along(assignment@transAssign))
            
            if (isTRUE(assignment@log$verbose)) {
              t_end_time <- proc.time()
              elapsed <- elapsedTime(t_start_time,
                                     t_end_time)
              message(blue('Transformation assignment '),'\t',green(cli::symbol$tick),' ',elapsed)
            }
            
            return(assignment)
          }
)
