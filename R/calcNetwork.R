#' @importFrom dplyr group_by summarise group_indices

calcNetwork <- function(MF,rel) {
  MF <- MF %>%
    tbl_df() %>%
    select(RetentionTime:Score)
  Nodes <- suppressWarnings( MF 
                             %>%
                              group_by(MF) %>%
                              summarise(Nodes = n()))
  Edges <- rel %>%
    group_by(MF) %>%
    summarise(Edges = n(),EdgeWeight = mean(r))
  
  MFs <- inner_join(Nodes,Edges, by = c('MF' = 'MF')) %>%
    mutate(Degree = Nodes * Edges)
  
  MF <- inner_join(MF,MFs, by = c('MF' = 'MF')) %>% 
  arrange(desc(Degree),desc(Edges),MF) %>%
    mutate(Cluster = NA)
  
  for (i in 1:length(unique(MF$MF))) {
    mf <- MF %>% 
      filter(MF == unique(MF)[i])
    clus <- filter(MF,`Measured m/z` %in% mf$`Measured m/z` & 
                     Isotope %in% mf$Isotope & 
                     Adduct %in% mf$Adduct &
                     `Measured M` %in% mf$`Measured M` &
                     Nodes %in% mf$Nodes &
                     Edges %in% mf$Edges &
                     Degree %in% mf$Degree)
    if (NA %in% clus$Cluster) {
      MF <- mutate(MF,Cluster = ifelse(`Measured m/z` %in% clus$`Measured m/z` & 
                                         Isotope %in% clus$Isotope & 
                                         Adduct %in% clus$Adduct &
                                         `Measured M` %in% clus$`Measured M` &
                                         Nodes %in% clus$Nodes &
                                         Edges %in% clus$Edges &
                                         Degree %in% clus$Degree,i,Cluster))
    }
  }
  cl <- group_by(MF,Cluster) %>% group_indices()
  MF <- mutate(MF,Cluster = cl) %>%
    arrange(Cluster)
  
  addIsoScores <- MF %>%
    group_by(Cluster,MF) %>% 
    summarise(AddIsoScore = addIsoScore(Adduct,Isotope,addRank = parameters@adducts,isoRank = parameters@isotopes))
  
  MF <- inner_join(MF,addIsoScores,by = c('Cluster','MF')) %>%
    mutate(`Average AddIsoScore` = AddIsoScore / Nodes,
           Plausibility = `Average AddIsoScore` / Edges)
  
  return(MF)
}