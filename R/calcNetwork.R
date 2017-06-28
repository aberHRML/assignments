
calcNetwork <- function(MF,rel) {
  MF <- tbl_df(select(MF,MF:Score))
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
  return(MF)
}