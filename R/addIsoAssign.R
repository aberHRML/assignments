
setGeneric("addIsoAssign", function(assignment)
  standardGeneric("addIsoAssign")
)

#' @importFrom dplyr arrange rowwise slice_sample left_join ungroup
#' @importFrom stringr str_detect
#' @importFrom mzAnnotation calcM ipMF
#' @importFrom igraph vertex.attributes V
#' @importFrom furrr furrr_options
#' @importFrom methods as

setMethod('addIsoAssign',signature = 'Assignment',
          function(assignment){
            
            if (isTRUE(assignment@log$verbose)) {
              startTime <- proc.time()
              message(blue('Adduct & isotopic assignment '),
                      cli::symbol$continue,'\r',
                      appendLF = FALSE)
            }
            
            assignment_technique <- technique(assignment)
            
            rel <- assignment %>% 
              relationships()
            
            if (str_detect(assignment_technique,'LC')){
              rel <- rel %>% 
                filter(RetentionTimeDiff <= assignment@RT_diff_limit)
            }
            
            rel <- rel %>% 
              filter(is.na(Transformation1) & 
                       is.na(Transformation2) & 
                       coefficient > 0) %>%
              filter(!(is.na(Isotope1) & 
                         !is.na(Isotope2) & 
                         Adduct1 == Adduct2 & 
                         log2IntensityRatio < 0)) %>%
              filter(!(!is.na(Isotope1) & 
                         is.na(Isotope2) & 
                         Adduct1 == Adduct2)) 
            
            M <- collateM(rel,
                          maxM(assignment))
            
            MFs <- generateMFs(M,
                               ppm(assignment),
                               MFrankThreshold(assignment),
                               adductRules(assignment),
                               isotopeRules(assignment),
                               AIS(assignment))
            
            graph_edges <- rel %>% 
              addMFs(MFs) %>%
              filter(MF1 == MF2) %>%
              mutate(RetentionTime1 = as.numeric(RetentionTime1),
                     RetentionTime2 = as.numeric(RetentionTime2)) %>%
              addNames()
            
            graph_nodes <- collateMFs(graph_edges,MFs)
            
            graph <- calcComponents(graph_nodes,
                                    graph_edges,
                                    assignment)
            
            counter <- 0
            
            assignment@addIsoAssign <- list()
            
            repeat {
              
              counter <- counter + 1
              
              message(paste0('iteration ',counter))
              
              if (counter > 1){
                graph <- assignment@addIsoAssign[[counter - 1]]$graph %>% 
                  activate(nodes) %>% 
                  dplyr::anti_join(assignment %>% 
                                     assignments() %>% 
                                     select(dplyr::any_of(c('Feature','Isotope','Adduct','MF'))),
                                   by = 'Feature') 
                
                if (length(graph) == 0) break()
                
                graph <- graph %>% 
                  clean(adductRules(assignment)) %>% 
                  recalcComponents(assignment)
              }
              
              filtered_graph <- graph 
              
              filtered_graph <- filtered_graph %>% 
                filterComponents(assignment,
                                 filters = componentFilters())
              
              assigned_features <- filtered_graph %>% 
                nodes() %>% 
                rename(Name = name) %>%
                mutate(Mode = str_sub(Feature,1,1))
              
              if (nrow(assigned_features) == 0){
                break()
              }
              
              assignment@addIsoAssign[[counter]] <- list(
                graph = graph,
                filtered_graph = filtered_graph,
                assigned = assigned_features
              )
              
              assignment@assignments <- bind_rows(
                assignment@assignments,
                assigned_features %>%
                  select(Name:`MF Plausibility (%)`,Mode) %>%
                  mutate(Iteration = paste0('A&I',counter))
              )
              
            }
            
            if (assignment@log$verbose == TRUE) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              message(blue('Adduct & isotopic assignment '),'\t',green(cli::symbol$tick),' ',elapsed)
            }
            
            return(assignment)
          })
