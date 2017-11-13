#' @importFrom dplyr arrange rowwise sample_n left_join
#' @importFrom mzAnnotation calcM calcMZ ppmError

setMethod('addIsoAssign',signature = 'Assignment',
          function(assignment){
            parameters <- assignment@parameters
            rel <- assignment@relationships %>% 
              filter(is.na(Transformation1) & is.na(Transformation2) & r > 0) %>%
              filter(!(is.na(Isotope1) & !is.na(Isotope2) & Adduct1 == Adduct2 & log2IntensityRatio < 0)) %>%
              filter(!(!is.na(Isotope1) & is.na(Isotope2) & Adduct1 == Adduct2)) 
            
            M <- bind_rows(select(rel,mz = `m/z1`,rt = rt1,Isotope = Isotope1, Adduct = Adduct1),
                           select(rel,mz = `m/z2`,rt = rt2,Isotope = Isotope2, Adduct = Adduct2)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz) %>%
              rowwise() %>%
              mutate(M = calcM(mz,adduct = Adduct,isotope = Isotope)) %>% 
              arrange(M) %>%
              filter(M <= parameters@maxM)
            
            clus <- makeCluster(parameters@nCores)
            MF <- sample_n(M,nrow(M)) %>%
              rowwise() %>% 
              parApply(cl = clus,1,function(x){MFgen(as.numeric(x[5]),as.numeric(x[1]),ppm = parameters@ppm)}) %>% 
              bind_rows() %>%
              left_join(M,by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
              rowwise() %>%
              mutate(`Theoretical m/z` = calcMZ(`Theoretical M`,Adduct,Isotope), 
                     `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)) %>%
              select(rt,MF,Isotope,Adduct,`Theoretical M`,`Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM Error`) %>%
              rowwise() %>%
              mutate(Score = MFscore(MF)) %>%
              filter(Score <= parameters@maxMFscore)
            stopCluster(clus)
            
            rel <- rel %>% addMFs(MF) %>%
              mutate(rt1 = as.numeric(rt1),rt2 = as.numeric(rt2))
            
            MFs <- bind_rows(select(rel,mz = `m/z1`,rt = rt1,Isotope = Isotope1, Adduct = Adduct1, MF = MF),
                             select(rel,mz = `m/z2`,rt = rt2,Isotope = Isotope2, Adduct = Adduct2,MF = MF)) %>%
              filter(!duplicated(.)) %>%
              mutate(rt = as.numeric(rt)) %>%
              arrange(mz)
            
            MF <- semi_join(MF,MFs,by = c('rt' = 'rt','MF' = 'MF','Isotope' = 'Isotope','Adduct' = 'Adduct','Measured m/z' = 'mz'))
            
            MF <- calcNetwork(MF,rel)
            
            filteredMF <- group_by(MF,Cluster) %>% 
              filter(Score == min(Score)) %>%
              group_by(`Measured m/z`) %>% 
              filter(Connectivity == max(Connectivity))
            
            filteredRel <- semi_join(rel,filteredMF,by = c('MF' = 'MF','Isotope1' = 'Isotope','Isotope2' = 'Isotope','Adduct1' = 'Adduct','Adduct2' = 'Adduct','m/z1' = 'Measured m/z','m/z2' = 'Measured m/z','rt1' = 'rt', 'rt2' = 'rt'))
            
            filteredMF <- calcNetwork(filteredMF,filteredRel) %>% 
              filter(Connectivity > 1) %>%
              group_by(`Measured m/z`) %>% 
              filter(Connectivity == max(Connectivity)) 
            
            addIsoScores <- filteredMF %>%
              group_by(Cluster) %>% 
              summarise(AddIsoScore = addIsoScore(Adduct,Isotope,addRank = parameters@adducts,isoRank = parameters@isotopes))
            
            filteredMF <- inner_join(filteredMF,addIsoScores,by = c('Cluster' = 'Cluster')) 
            
            filteredMF <- filteredMF %>%
              group_by(`Measured m/z`) %>%
              filter(AddIsoScore == min(AddIsoScore))
            
            filteredRel <- semi_join(filteredRel,filteredMF,by = c('MF' = 'MF','Isotope1' = 'Isotope','Isotope2' = 'Isotope','Adduct1' = 'Adduct','Adduct2' = 'Adduct','m/z1' = 'Measured m/z','m/z2' = 'Measured m/z','rt1' = 'rt', 'rt2' = 'rt'))
            
            filteredMF <- calcNetwork(filteredMF,filteredRel) %>%
              filter(Connectivity > 1) 
            
            addIsoScores <- filteredMF %>%
              group_by(Cluster) %>% 
              summarise(AddIsoScore = addIsoScore(Adduct,Isotope,addRank = parameters@adducts,isoRank = parameters@isotopes))
            
            filteredMF <- inner_join(filteredMF,addIsoScores,by = c('Cluster' = 'Cluster')) 
            
            filteredMF <- filteredMF %>%
              group_by(`Measured m/z`) %>% 
              filter(Score == min(Score)) %>%
              mutate(absPPM = abs(`PPM Error`)) %>%
              filter(absPPM == min(absPPM)) %>% 
              select(rt:AddIsoScore)
            
            filteredMF <- calcNetwork(filteredMF,filteredRel) %>%
              filter(Connectivity > 1)
            
            adducts <- lapply(parameters@adducts,function(y){tibble(Adduct = y)})
            adducts <- bind_rows(adducts,.id = 'Mode')
            
            assigned <- select(filteredMF,rt:Score) %>%
              inner_join(adducts,c('Adduct' = 'Adduct')) %>%
              arrange(`MF`)
            
            assignment@assignments <- assigned
            assignment@addIsoAssign <- list(MFs = MF, relationships = rel, filteredMFs = filteredMF, filteredRelationships = filteredRel,assigned = assigned)
            return(assignment)
          })



