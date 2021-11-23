#' @importFrom stringr str_c

setGeneric("transformationAssign", function(assignment)
  standardGeneric("transformationAssign"))

#' @importFrom dplyr full_join select distinct
#' @importFrom mzAnnotation transformMF

setMethod('transformationAssign',signature = 'Assignment',
          function(assignment){
            
            parameters <- as(assignment,'AssignmentParameters')
            count <- length(assignment@transAssign)
            assigned <- assignment@assignments
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              message(blue(str_c('Transformation assignment iteration ', count + 1,' ')),cli::symbol$continue,'\r',appendLF = FALSE)
            }
            
            rel <- assignment@relationships %>%
              filter((`m/z1` %in% assigned$`Measured m/z` | (`m/z2` %in% assigned$`Measured m/z`)) & !(`m/z1` %in% assigned$`Measured m/z` & (`m/z2` %in% assigned$`Measured m/z`)))
            
            mz1 <- rel %>%
              semi_join(assigned,by = c('m/z1' = 'Measured m/z', 'Adduct1' = 'Adduct', 'Isotope1' = 'Isotope')) %>%
              filter(is.na(Transformation1))
            mz2 <- rel %>%
              semi_join(assigned,by = c('m/z2' = 'Measured m/z', 'Adduct2' = 'Adduct', 'Isotope2' = 'Isotope')) %>%
              filter(is.na(Transformation2))
            
            rel <- bind_rows(mz1,mz2)
            
            if (nrow(rel) > 0) {
              M <- bind_rows(select(rel,mz = `m/z1`,RetentionTime = RetentionTime1,Isotope = Isotope1, Adduct = Adduct1, Feature = Feature1),
                             select(rel,mz = `m/z2`,RetentionTime = RetentionTime2,Isotope = Isotope2, Adduct = Adduct2, Feature = Feature2)) %>%
                distinct() %>%
                arrange(mz) %>%
                rowwise() %>%
                mutate(M = calcM(mz,Adduct,Isotope)) %>% 
                arrange(M) %>%
                filter(M <= parameters@maxM) %>%
                filter(!(mz %in% assigned$`Measured m/z`))
              
              nM <- nrow(M)
              
              MF <- M %>%
                ungroup() %>%
                slice_sample(n = nM) %>%
                split(1:nrow(.)) %>%
                future_map(~{
                  mf <- MFgen(.x$M,.x$mz,ppm = parameters@ppm) 
                  
                  if (nrow(mf) > 0) {
                    mf %>%
                      left_join(M,by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
                      rowwise() %>%
                      mutate(`Theoretical m/z` = calcMZ(`Theoretical M`,Adduct,Isotope), 
                             `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)) %>%
                      select(Feature,RetentionTime,MF,Isotope,Adduct,`Theoretical M`,
                             `Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM Error`) %>%
                      rowwise() %>%
                      mutate(Score = MFscore(MF),
                             `PPM Error` = abs(`PPM Error`),
                             AddIsoScore = addIsoScore(Adduct,
                                                       Isotope,
                                                       parameters@adducts,
                                                       parameters@isotopes)) %>%
                      ungroup() %>%
                      filter(Score == min(Score,na.rm = TRUE)) %>%
                      filter(Score < parameters@maxMFscore)
                  } else {
                    return(NULL)
                  }
                },.options = furrr_options(seed = 1234)) %>% 
                bind_rows()
              
              if (nrow(MF) > 0) {
                
                MF <- MF %>%
                  bind_rows(assigned %>%
                              select(names(MF)[!(names(MF) == 'AddIsoScore')]) %>%
                              rowwise() %>%
                              mutate(AddIsoScore = addIsoScore(Adduct,Isotope,parameters@adducts,parameters@isotopes)))
                rel <- rel %>% 
                  addMFs(MF,identMF = F) %>%
                  mutate(RetentionTime1 = as.numeric(RetentionTime1),RetentionTime2 = as.numeric(RetentionTime2)) %>%
                  addNames()
                
                if (nrow(rel) > 0) {
                  MFs <- bind_rows(select(rel,Name = Name1,Feature = Feature1,mz = `m/z1`,RetentionTime = RetentionTime1,Isotope = Isotope1, Adduct = Adduct1, MF = MF1),
                                   select(rel,Name = Name2,Feature = Feature2,mz = `m/z2`,RetentionTime = RetentionTime2,Isotope = Isotope2, Adduct = Adduct2,MF = MF2)) %>%
                    mutate(RetentionTime = as.numeric(RetentionTime)) %>%
                    arrange(mz) %>%
                    select(-mz) %>%
                    left_join(MF, by = c("Feature", "RetentionTime", "Isotope", "Adduct",'MF')) %>%
                    distinct() %>%
                    mutate(ID = 1:nrow(.))
                  
                  graph <- calcComponents(MFs,rel,parameters)
                  
                  filters <- tibble(Measure = c('Plausibility','Size','AIS','Score','PPM Error'),
                                    Direction = c(rep('max',3),rep('min',2)))
                  
                  filteredGraph <- graph
                  
                  for (i in 1:nrow(filters)) { 
                    f <- filters[i,]
                    filteredGraph <- filteredGraph %>%
                      activate(nodes) %>%
                      filter(name %in% {filteredGraph %>% 
                          vertex.attributes() %>% 
                          as_tibble() %>%
                          eliminate(f$Measure,f$Direction) %>%
                          .$name}) 
                    if (V(filteredGraph) %>% length() > 0) {
                      filteredGraph <- filteredGraph %>%
                        recalcComponents(parameters)
                    } else {
                      break()
                    }
                  }
                  
                  newlyAssigned <- filteredGraph %>%
                    vertex.attributes() %>% 
                    as_tibble() %>%
                    rename(Name = name) %>%
                    mutate(Mode = str_sub(Feature,1,1)) %>%
                    filter(!(Name %in% assigned$Name)) %>%
                    select(Name:Score,Mode) %>%
                    mutate(Iteration = str_c('T',count + 1))
                  
                  outputs <- list(
                    graph = graph,
                    filteredGraph = filteredGraph,
                    assigned = newlyAssigned)
                  
                  assignment@assignments <- bind_rows(assignment@assignments,newlyAssigned)
                  
                  if (count == 0) {
                    assignment@transAssign <- list(`1` = outputs)
                  } else {
                    assignment@transAssign <- c(assignment@transAssign,list(outputs))
                  }
                } else {
                  assignment@transAssign <- c(assignment@transAssign,list(list()))
                }
              } else {
                assignment@transAssign <- c(assignment@transAssign,list(list())) 
              }
            } else {
              assignment@transAssign <- c(assignment@transAssign,list(list()))
            }
            names(assignment@transAssign)[count + 1] <- count + 1
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              message(blue(str_c('Transformation assignment iteration ', count + 1,' ')),'\t',green(cli::symbol$tick),' ',elapsed)
            }
            
            return(assignment)
          }
)
