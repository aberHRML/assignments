
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
              message(blue('Adduct & isotopic assignment '),cli::symbol$continue,'\r',appendLF = FALSE)
            }
            
            rel <- assignment %>% 
              relationships() %>% 
              filter(is.na(Transformation1) & 
                       is.na(Transformation2) & 
                       r > 0) %>%
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
                              adducts(assignment),
                              adductRules(assignment),
                              isotopes(assignment),
                              isotopeRules(assignment))
              
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
  
            filtered_graph <- filterComponents(graph,
                                               assignment)
            
            assignment@addIsoAssign <- list(
              graph = graph,
              filtered_graph = filtered_graph,
              assigned = filtered_graph %>% 
                nodes() %>% 
                rename(Name = name) %>%
                mutate(Mode = str_sub(Feature,1,1))
            )
            
            assignment@assignments <- assignment@addIsoAssign$assigned %>%
              select(Name:`MF Plausibility (%)`,Mode) %>%
              mutate(Iteration = 'A&I')
            
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
