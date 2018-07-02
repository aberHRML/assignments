#' @importFrom stringr str_c
#' @importFrom dplyr full_join select
#' @importFrom mzAnnotation transformMF

setMethod('transformationAssign',signature = 'Assignment',
          function(assignment){
            parameters <- assignment@parameters
            count <- length(assignment@transAssign)
            assigned <- assignment@assignments
            
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
              M <- bind_rows(select(rel,mz = `m/z1`,RetentionTime = RetentionTime1,Isotope = Isotope1, Adduct = Adduct1),
                             select(rel,mz = `m/z2`,RetentionTime = RetentionTime2,Isotope = Isotope2, Adduct = Adduct2)) %>%
                filter(!duplicated(.)) %>%
                arrange(mz) %>%
                rowwise() %>%
                mutate(M = calcM(mz,Adduct,Isotope)) %>% 
                arrange(M) %>%
                filter(M <= parameters@maxM) %>%
                filter(!(mz %in% assigned$`Measured m/z`))
              
              clus <- makeCluster(parameters@nCores)
              MF <- sample_n(M,nrow(M)) %>%
                rowwise() %>% 
                parApply(cl = clus,1,function(x,ppm){MFassign:::MFgen(as.numeric(x[5]),as.numeric(x[1]),ppm)},ppm = parameters@ppm) %>% 
                bind_rows() %>%
                left_join(M,by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
                rowwise() %>%
                mutate(`Theoretical m/z` = calcMZ(`Theoretical M`,Adduct,Isotope), 
                       `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)) %>%
                select(RetentionTime,MF,Isotope,Adduct,`Theoretical M`,`Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM Error`) %>%
                rowwise() %>%
                mutate(Score = MFscore(MF)) %>%
                filter(Score <= parameters@maxMFscore)
              stopCluster(clus)
              
              MFs <- select(MF,RetentionTime,MF,Isotope,Adduct,`Measured m/z`) %>%
                bind_rows(select(assigned,RetentionTime,MF,Isotope,Adduct,`Measured m/z`))
              
              rel <- filter(rel, `m/z1` %in% MF$`Measured m/z` | `m/z2` %in% MF$`Measured m/z`) %>%
                left_join(MFs,by = c('RetentionTime1' = 'RetentionTime','m/z1' = 'Measured m/z')) %>%
                rename(MF1 = MF)
              rel[is.na(rel)] <- ''
              
              rel <- rel %>%
                mutate(RetentionTime1 = as.numeric(RetentionTime1),RetentionTime2 = as.numeric(RetentionTime2))
              
              rel <- filter(rel,Isotope1 == Isotope & Adduct1 == Adduct) %>%
                select(`m/z1`:MF1) %>%
                left_join(MFs,by = c('RetentionTime2' = 'RetentionTime','m/z2' = 'Measured m/z')) %>%
                rename(MF2 = MF)
              rel[is.na(rel)] <- ''
              
              rel <- filter(rel,Isotope2 == Isotope & Adduct2 == Adduct) %>%
                select(-(Isotope:Adduct))
              rel[rel == ''] <- NA
              
              rel <- rel %>%
                mutate(RetentionTime1 = as.numeric(RetentionTime1),RetentionTime2 = as.numeric(RetentionTime2))
              
              trans1 <- rel %>%
                select(MF2,Transformation1) %>%
                rowwise() %>%
                mutate(TransformedMF = transformMF(MF2,Transformation1))
              
              trans2 <- rel %>%
                select(MF1,Transformation2) %>%
                rowwise() %>%
                mutate(TransformedMF = transformMF(MF1,Transformation2))
              
              rel <- rel %>%
                mutate(TransformedMF1 = trans1$TransformedMF, TransformedMF2 = trans2$TransformedMF) %>%
                filter(TransformedMF1 == TransformedMF2)
              
              MFs <- bind_rows(select(rel,RetentionTime = RetentionTime1,mz = `m/z1`,Isotope = Isotope1, Adduct = Adduct1, MF = MF1),
                               select(rel,RetentionTime = RetentionTime2,mz = `m/z2`,Isotope = Isotope2, Adduct = Adduct2,MF = MF2)) %>%
                filter(!duplicated(.)) %>%
                arrange(mz) %>%
                filter(!(mz %in% assigned$`Measured m/z`))
              
              MF <- semi_join(MF,MFs,by = c('RetentionTime' = 'RetentionTime','MF' = 'MF','Isotope' = 'Isotope','Adduct' = 'Adduct','Measured m/z' = 'mz'))
              
              MF <- calcNetworktrans(MF,rel)
              
              filteredMF <- group_by(MF,Cluster) %>% 
                filter(Score == min(Score)) %>%
                group_by(`Measured m/z`) %>% 
                filter(Degree == max(Degree))
              
              filteredRel1 <- semi_join(rel,filteredMF,by = c('RetentionTime1' = 'RetentionTime','MF1' = 'MF','Isotope1' = 'Isotope','Adduct1' = 'Adduct','m/z1' = 'Measured m/z'))
              filteredRel2 <- semi_join(rel,filteredMF,by = c('RetentionTime2' = 'RetentionTime','MF2' = 'MF','Isotope2' = 'Isotope','Adduct2' = 'Adduct','m/z2' = 'Measured m/z'))
              filteredRel <- bind_rows(filteredRel1,filteredRel2) %>%
                select(`m/z1`:MF2)
              
              filteredMF <- calcNetworktrans(filteredMF,filteredRel) %>%
                group_by(`Measured m/z`) %>% 
                filter(Degree == max(Degree))
              
              addIsoScores <- filteredMF %>%
                group_by(Cluster) %>% 
                summarise(AddIsoScore = addIsoScore(Adduct,Isotope,addRank = parameters@adducts,isoRank = parameters@isotopes))
              
              filteredMF <- inner_join(filteredMF,addIsoScores,by = c('Cluster' = 'Cluster')) 
              
              filteredMF <- filteredMF %>%
                group_by(`Measured m/z`) %>%
                filter(AddIsoScore == min(AddIsoScore))
              
              filteredRel1 <- semi_join(rel,filteredMF,by = c('RetentionTime1' = 'RetentionTime','MF1' = 'MF','Isotope1' = 'Isotope','Adduct1' = 'Adduct','m/z1' = 'Measured m/z'))
              filteredRel2 <- semi_join(rel,filteredMF,by = c('RetentionTime2' = 'RetentionTime','MF2' = 'MF','Isotope2' = 'Isotope','Adduct2' = 'Adduct','m/z2' = 'Measured m/z'))
              filteredRel <- bind_rows(filteredRel1,filteredRel2) 
              
              filteredMF <- calcNetworktrans(filteredMF,filteredRel)
              
              addIsoScores <- filteredMF %>%
                group_by(Cluster) %>% 
                summarise(AddIsoScore = addIsoScore(Adduct,Isotope,addRank = parameters@adducts,isoRank = parameters@isotopes))
              
              filteredMF <- inner_join(filteredMF,addIsoScores,by = c('Cluster' = 'Cluster')) 
              
              filteredMF <- filteredMF %>%
                group_by(`Measured m/z`) %>% 
                filter(Score == min(Score)) %>%
                mutate(absPPM = abs(`PPM Error`)) %>%
                filter(absPPM == min(absPPM)) %>% 
                select(RetentionTime:AddIsoScore)
              
              filteredMF <- calcNetworktrans(filteredMF,filteredRel)
              
              adducts <- lapply(parameters@adducts,function(y){tibble(Adduct = y)})
              adducts <- bind_rows(adducts,.id = 'Mode')
              
              assignments <- select(filteredMF,RetentionTime:Score) %>%
                inner_join(adducts,c('Adduct' = 'Adduct')) %>%
                arrange(`MF`)
              assignment@assignments <- bind_rows(assignment@assignments,assignments)
              
              if (count == 0) {
                assignment@transAssign <- list(`1` = list(MFs = MF, relationships = rel, filteredMFs = filteredMF, filteredRelationships = filteredRel,assigned = assigned))
              } else {
                assignment@transAssign <- c(assignment@transAssign,list(list(MFs = MF, relationships = rel, filteredMFs = filteredMF, filteredRelationships = filteredRel,assigned = assigned)))
              }  
            } else {
              assignment@transAssign <- c(assignment@transAssign,list(list()))
              
            }
            names(assignment@transAssign)[count + 1] <- count + 1
            return(assignment)
          }
)
