#' @importFrom stringr str_c
#' @importFrom dplyr full_join select distinct
#' @importFrom mzAnnotation transformMF

setMethod('transformationAssign',signature = 'Assignment',
          function(assignment){
            
            parameters <- assignment@parameters
            count <- length(assignment@transAssign)
            assigned <- assignment@assignments
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              cat(blue(str_c('Transformation assignment iteration ', count + 1,' ')),cli::symbol$continue,'\r',sep = '')
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
              
              slaves <- {nrow(M) / 20} %>%
                round()
              
              if (slaves > parameters@nCores) {
                slaves <- parameters@nCores
              }
                
              clus <- makeCluster(slaves)
              
              MF <- sample_n(M,nM) %>%
                split(1:nrow(.)) %>%
                parLapply(cl = clus,function(x,parameters){MFgen(x$M,x$mz,ppm = parameters@ppm)},parameters = parameters) %>% 
                bind_rows() %>%
                left_join(M,by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
                rowwise() %>%
                mutate(`Theoretical m/z` = calcMZ(`Theoretical M`,Adduct,Isotope), 
                       `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)) %>%
                select(Feature,RetentionTime,MF,Isotope,Adduct,`Theoretical M`,`Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM Error`) %>%
                rowwise() %>%
                mutate(Score = MFscore(MF)) %>%
                filter(Score <= parameters@maxMFscore)
              stopCluster(clus)
             
              if (nrow(MF) > 0) {
                MF <- MF %>%
                  bind_rows(assigned %>%
                              select(names(MF)))
                rel <- rel %>% 
                  addMFs(MF,identMF = F) %>%
                  mutate(RetentionTime1 = as.numeric(RetentionTime1),RetentionTime2 = as.numeric(RetentionTime2)) %>%
                  addNames()
                
                MFs <- bind_rows(select(rel,Name = Name1,Feature = Feature1,mz = `m/z1`,RetentionTime = RetentionTime1,Isotope = Isotope1, Adduct = Adduct1, MF = MF1),
                                 select(rel,Name = Name2,Feature = Feature2,mz = `m/z2`,RetentionTime = RetentionTime2,Isotope = Isotope2, Adduct = Adduct2,MF = MF2)) %>%
                  distinct() %>%
                  mutate(RetentionTime = as.numeric(RetentionTime)) %>%
                  arrange(mz) %>%
                  select(-mz) %>%
                  left_join(MF, by = c("Feature", "RetentionTime", "Isotope", "Adduct",'MF')) %>%
                  mutate(ID = 1:nrow(.)) %>%
                  rowwise() %>%
                  mutate(AddIsoScore = addIsoScore(Adduct,Isotope,parameters@adducts,parameters@isotopes),
                         `PPM Error` = abs(`PPM Error`)) %>%
                  tbl_df()
                
                graph <- calcComponents(MFs,rel)
                
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
                      recalcComponents()
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
            names(assignment@transAssign)[count + 1] <- count + 1
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              cat(blue(str_c('Transformation assignment iteration ', count + 1,' ')),'\t',green(cli::symbol$tick),' ',elapsed,'\n',sep = '')
            }
            
            return(assignment)
          }
)
