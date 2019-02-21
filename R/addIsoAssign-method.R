#' @importFrom dplyr arrange rowwise sample_n left_join
#' @importFrom stringr str_detect
#' @importFrom mzAnnotation calcM calcMZ ppmError
#' @importFrom igraph vertex.attributes
#' @importFrom parallel parLapply

setMethod('addIsoAssign',signature = 'Assignment',
          function(assignment){
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              cat(blue('Adduct & isotope assignment '),cli::symbol$continue,'\r',sep = '')
            }
            
            parameters <- assignment@parameters
            
            rel <- assignment@relationships %>% 
              filter(is.na(Transformation1) & is.na(Transformation2) & r > 0) %>%
              filter(!(is.na(Isotope1) & !is.na(Isotope2) & Adduct1 == Adduct2 & log2IntensityRatio < 0)) %>%
              filter(!(!is.na(Isotope1) & is.na(Isotope2) & Adduct1 == Adduct2)) 
            
            if (str_detect(parameters@technique,'-LC')) {
              rel <- rel %>%
                mutate(rtDiff = abs(RetentionTime1 - RetentionTime2)) %>%
                filter(rtDiff <= parameters@RTwindow) %>%
                select(-rtDiff)
            }
            
            M <- bind_rows(select(rel,mz = `m/z1`,RetentionTime = RetentionTime1,Isotope = Isotope1, Adduct = Adduct1, Feature = Feature1),
                           select(rel,mz = `m/z2`,RetentionTime = RetentionTime2,Isotope = Isotope2, Adduct = Adduct2, Feature = Feature2)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz) %>%
              rowwise() %>%
              mutate(M = calcM(mz,adduct = Adduct,isotope = Isotope)) %>% 
              arrange(M) %>%
              filter(M <= parameters@maxM)
            
            nM <- nrow(M)
            
            clus <- makeCluster(parameters@nCores)

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
            
            rel <- rel %>% 
              addMFs(MF) %>%
              filter(MF1 == MF2) %>%
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
                     `PPM Error` = abs(`PPM Error`))
            
            graph <- calcComponents(MFs,rel)
            
            filters <- tibble(Measure = c('Plausibility','Size','AddIsoScore','Score','PPM Error'),
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
                    .$name}) %>%
                recalcComponents()
            }
            
            assignment@addIsoAssign <- list(
              graph = graph,
              filteredGraph = filteredGraph,
              assigned = filteredGraph %>% 
                vertex.attributes() %>% 
                as_tibble()
            )
            
            assignment@assignments <- assignment@addIsoAssign$assigned
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              cat(blue('Adduct & isotope assignment '),'\t',green(cli::symbol$tick),' ',elapsed,'\n',sep = '')
            }
            
            return(assignment)
          })
