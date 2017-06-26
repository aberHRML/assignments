
setMethod('addIsoAssign',signature = 'Annotation',
          function(x){
            parameters <- x@parameters
            rel <- x@relationships %>% 
              filter(is.na(Transformation1) & is.na(Transformation2))
            
            M <- bind_rows(select(rel,mz = mz1,Isotope = Isotope1, Adduct = Adduct1),
                        select(rel,mz = mz2,Isotope = Isotope2, Adduct = Adduct2)) %>%
              filter(!duplicated(.)) %>%
              arrange(mz) %>%
              rowwise() %>%
              mutate(M = calcM(mz,Isotope,Adduct)) %>% 
              arrange(M)
            
            MF <- rowwise(M) %>% 
              apply(1,function(x){MFgen(as.numeric(x[4]),as.numeric(x[1]),ppm = parameters@ppm)}) %>% 
              bind_rows() %>%
              tbl_df() %>% 
              left_join(M,by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
              rowwise() %>%
              mutate(`Theoretical m/z` = calcMZ(`Theoretical M`,Isotope,Adduct), `PPM Error` = round((`Measured m/z` - `Theoretical m/z`)/`Theoretical m/z` * 10^6,5)) %>%
              select(MF,Isotope,Adduct,`Theoretical M`,`Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM Error`) %>%
              rowwise() %>%
              mutate(Score = MFscore(`MF`))
            
            filteredRel <- filter(rel, mz1 %in% MF$`Measured m/z`,mz2 %in% MF$`Measured m/z`)
          })