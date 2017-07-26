#' @importFrom dplyr arrange rowwise sample_n left_join

setMethod('addIsoAssign',signature = 'Assignment',
          function(x){
            parameters <- x@parameters
            rel <- x@relationships %>% 
              filter(is.na(Transformation1) & is.na(Transformation2) & r > 0) %>%
              filter(!(is.na(Isotope1) & !is.na(Isotope2) & Adduct1 == Adduct2 & log2IntensityRatio < 0)) %>%
              filter(!(!is.na(Isotope1) & is.na(Isotope2) & Adduct1 == Adduct2)) 
            
            M <- bind_rows(select(rel,mz = mz1,Isotope = Isotope1, Adduct = Adduct1),
                           select(rel,mz = mz2,Isotope = Isotope2, Adduct = Adduct2)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz) %>%
              rowwise() %>%
              mutate(M = calcM(mz,Isotope,Adduct)) %>% 
              arrange(M) %>%
              filter(M <= parameters@maxM)
            
            clus <- makeCluster(parameters@nCores)
            MF <- sample_n(M,nrow(M)) %>%
              rowwise() %>% 
              parApply(cl = clus,1,function(x){MFgen(as.numeric(x[4]),as.numeric(x[1]),ppm = parameters@ppm)}) %>% 
              bind_rows() %>%
              tbl_df() %>% 
              left_join(M,by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
              rowwise() %>%
              mutate(`Theoretical m/z` = calcMZ(`Theoretical M`,Isotope,Adduct), `PPM Error` = round((`Measured m/z` - `Theoretical m/z`)/`Theoretical m/z` * 10^6,5)) %>%
              select(MF,Isotope,Adduct,`Theoretical M`,`Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM Error`) %>%
              rowwise() %>%
              mutate(Score = MFscore(MF)) %>%
              filter(Score <= parameters@maxMFscore)
            stopCluster(clus)
            
            rel <- addMFs(rel,MF)
            
            MFs <- bind_rows(select(rel,mz = mz1,Isotope = Isotope1, Adduct = Adduct1, MF = MF),
                             select(rel,mz = mz2,Isotope = Isotope2, Adduct = Adduct2,MF = MF)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz)
            
            MF <- semi_join(MF,MFs,by = c('MF' = 'MF','Isotope' = 'Isotope','Adduct' = 'Adduct','Measured m/z' = 'mz'))
            
           MF <- calcNetwork(MF,rel)
            
            filteredMF <- group_by(MF,Cluster) %>% 
              filter(Score == min(Score)) %>%
              group_by(`Measured m/z`) %>% 
              filter(Connectivity == max(Connectivity))
            
            filteredRel <- semi_join(rel,filteredMF,by = c('MF' = 'MF','Isotope1' = 'Isotope','Isotope2' = 'Isotope','Adduct1' = 'Adduct','Adduct2' = 'Adduct','mz1' = 'Measured m/z','mz2' = 'Measured m/z'))
            
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
              filter(AddIsoScore == max(AddIsoScore))
            
            filteredRel <- semi_join(filteredRel,filteredMF,by = c('MF' = 'MF','Isotope1' = 'Isotope','Isotope2' = 'Isotope','Adduct1' = 'Adduct','Adduct2' = 'Adduct','mz1' = 'Measured m/z','mz2' = 'Measured m/z'))
            
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
              select(MF:AddIsoScore)
            
            filteredMF <- calcNetwork(filteredMF,filteredRel) %>%
              filter(Connectivity > 1)
            
            adducts <- lapply(parameters@adducts,function(y){tibble(Adduct = y)})
            adducts <- bind_rows(adducts,.id = 'Mode')
            
            assigned <- select(filteredMF,MF:Score) %>%
              inner_join(adducts,c('Adduct' = 'Adduct')) %>%
              arrange(`MF`)
            
            x@assignments <- assigned
            x@addIsoAssign <- list(MFs = MF, relationships = rel, filteredMFs = filteredMF, filteredRelationships = filteredRel)
            return(x)
          })



