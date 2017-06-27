
setMethod('addIsoAssign',signature = 'Annotation',
          function(x){
            parameters <- x@parameters
            rel <- x@relationships %>% 
              filter(is.na(Transformation1) & is.na(Transformation2) & Isotope1 %in% c(NA,parameters@isotopes) & Isotope2 %in% c(NA,parameters@isotopes) & r > 0)
            
            M <- bind_rows(select(rel,mz = mz1,Isotope = Isotope1, Adduct = Adduct1),
                           select(rel,mz = mz2,Isotope = Isotope2, Adduct = Adduct2)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz) %>%
              rowwise() %>%
              mutate(M = calcM(mz,Isotope,Adduct)) %>% 
              arrange(M)
            
            clus <- makeCluster(parameters@nCores)
            clusterExport(clus,c('MFgen','generateMF','mutate','parameters'))
            MF <- rowwise(M) %>% 
              parApply(cl = clus,1,function(x){MFgen(as.numeric(x[4]),as.numeric(x[1]),ppm = parameters@ppm)}) %>% 
              bind_rows() %>%
              tbl_df() %>% 
              left_join(M,by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
              rowwise() %>%
              mutate(`Theoretical m/z` = calcMZ(`Theoretical M`,Isotope,Adduct), `PPM Error` = round((`Measured m/z` - `Theoretical m/z`)/`Theoretical m/z` * 10^6,5)) %>%
              select(MF,Isotope,Adduct,`Theoretical M`,`Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM Error`) %>%
              rowwise() %>%
              mutate(Score = MFscore(MF))
            stopCluster(clus)
            
            MFs <- select(MF,MF,Isotope,Adduct,`Measured m/z`,)
            
            rel <- filter(rel, mz1 %in% MF$`Measured m/z`,mz2 %in% MF$`Measured m/z`) %>%
              left_join(MFs,by = c('mz1' = 'Measured m/z')) %>%
              rename(MF1 = MF)
            rel[is.na(rel)] <- ''
            rel <- filter(rel,Isotope1 == Isotope & Adduct1 == Adduct) %>%
              select(mz1:MF1) %>%
              left_join(MFs,by = c('mz2' = 'Measured m/z')) %>%
              rename(MF2 = MF)
            rel[is.na(rel)] <- ''
            rel <- filter(rel,Isotope2 == Isotope & Adduct2 == Adduct) %>%
              filter(MF1 == MF2) %>%
              select(mz1:Isotope2,Adduct1:MF1) %>%
              rename(MF = MF1)
            rel[rel == ''] <- NA
            
            MFs <- bind_rows(select(rel,mz = mz1,Isotope = Isotope1, Adduct = Adduct1, MF = MF),
                             select(rel,mz = mz2,Isotope = Isotope2, Adduct = Adduct2,MF = MF)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz)
            
            MF <- semi_join(MF,MFs,by = c('MF' = 'MF','Isotope' = 'Isotope','Adduct' = 'Adduct','Measured m/z' = 'mz'))
            
            Nodes <- suppressWarnings(group_by(MF,MF)) %>%
              summarise(Nodes = n())
            Edges <- group_by(rel,MF) %>%
              summarise(Edges = n(),EdgeWeight = mean(r))
            
            MFs <- inner_join(Nodes,Edges, by = c('MF' = 'MF')) %>%
              mutate(Connectivity = Nodes * Edges)
            
            MF <- inner_join(MF,MFs, by = c('MF' = 'MF')) %>% 
              arrange(desc(Connectivity),desc(Edges),MF) %>%
              mutate(Cluster = NA)
            
            for (i in 1:length(unique(MF$MF))) {
              mf <- filter(MF,MF == unique(MF)[i])
              clus <- filter(MF,`Measured m/z` %in% mf$`Measured m/z` & 
                               Isotope %in% mf$Isotope & 
                               Adduct %in% mf$Adduct &
                               `Measured M` %in% mf$`Measured M` &
                               Nodes %in% mf$Nodes &
                               Edges %in% mf$Edges &
                               Connectivity %in% mf$Connectivity)
              if (NA %in% clus$Cluster) {
                MF <- mutate(MF,Cluster = ifelse(`Measured m/z` %in% clus$`Measured m/z` & 
                                                   Isotope %in% clus$Isotope & 
                                                   Adduct %in% clus$Adduct &
                                                   `Measured M` %in% clus$`Measured M` &
                                                   Nodes %in% clus$Nodes &
                                                   Edges %in% clus$Edges &
                                                   Connectivity %in% clus$Connectivity,i,Cluster))
              }
            }
            cl <- group_by(MF,Cluster) %>% group_indices()
            MF <- mutate(MF,Cluster = cl) %>%
              arrange(Cluster)
            
            filteredMF <- group_by(MF,Cluster) %>% filter(Score == max(Score))
            filteredMF <- group_by(MF,`Measured m/z`) %>% filter(Connectivity == max(Connectivity))
          })



