#' @importFrom dplyr arrange rowwise slice_sample left_join ungroup
#' @importFrom stringr str_detect
#' @importFrom mzAnnotation calcM calcMZ ppmError
#' @importFrom igraph vertex.attributes V
#' @importFrom furrr furrr_options
#' @importFrom methods as

setMethod('addIsoAssign',signature = 'Assignment',
          function(assignment){
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              message(blue('Adduct & isotope assignment '),cli::symbol$continue,'\r',appendLF = FALSE)
            }
            
            parameters <- as(assignment,'AssignmentParameters')
            
            rel <- assignment@relationships %>% 
              filter(is.na(Transformation1) & is.na(Transformation2) & r > 0) %>%
              filter(!(is.na(Isotope1) & !is.na(Isotope2) & Adduct1 == Adduct2 & log2IntensityRatio < 0)) %>%
              filter(!(!is.na(Isotope1) & is.na(Isotope2) & Adduct1 == Adduct2)) 
            
            M <- bind_rows(select(rel,mz = `m/z1`,RetentionTime = RetentionTime1,Isotope = Isotope1, Adduct = Adduct1, Feature = Feature1),
                           select(rel,mz = `m/z2`,RetentionTime = RetentionTime2,Isotope = Isotope2, Adduct = Adduct2, Feature = Feature2)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz) %>%
              rowwise() %>%
              mutate(M = calcM(mz,adduct = Adduct,isotope = Isotope)) %>% 
              arrange(M) %>%
              filter(M <= parameters@maxM)
            
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
              
            rel <- rel %>% 
              addMFs(MF) %>%
              filter(MF1 == MF2) %>%
              mutate(RetentionTime1 = as.numeric(RetentionTime1),RetentionTime2 = as.numeric(RetentionTime2)) %>%
              addNames()
            
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
            
            assignment@addIsoAssign <- list(
              graph = graph,
              filteredGraph = filteredGraph,
              assigned = filteredGraph %>% 
                vertex.attributes() %>% 
                as_tibble() %>%
                rename(Name = name) %>%
                mutate(Mode = str_sub(Feature,1,1))
            )
            
            assignment@assignments <- assignment@addIsoAssign$assigned %>%
              select(Name:Score,Mode) %>%
              mutate(Iteration = 'A&I')
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              message(blue('Adduct & isotope assignment '),'\t',green(cli::symbol$tick),' ',elapsed)
            }
            
            return(assignment)
          })
