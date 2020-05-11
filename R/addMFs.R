#' @importFrom dplyr rename

addMFs <- function(rel,MF,identMF = T){

  if (identMF == T) {
    relations <- rel %>%
      filter(Feature1 %in% MF$Feature, Feature2 %in% MF$Feature)  
  } else {
    relations <- rel %>%
      filter(Feature1 %in% MF$Feature | Feature2 %in% MF$Feature)  
  }
   relations <- relations %>%
    left_join(MF %>%
                select(Feature,MF,Isotope,Adduct,`Measured m/z`),by = c('Feature1' = 'Feature')) %>%
    rename(MF1 = MF)
   
   chr_columns <- relations %>%
     map_lgl(is.character)
   
   relations[,chr_columns] <- relations[,chr_columns] %>%
     {
       .[is.na(.)] <- ''
       return(.)
     }
  
  relations <- relations %>%
    filter(Isotope1 == Isotope & Adduct1 == Adduct) %>%
    select(-(Isotope:`Measured m/z`)) %>%
    left_join(MF %>%
                select(Feature,MF,Isotope,Adduct,`Measured m/z`),by = c('Feature2' = 'Feature')) %>%
    rename(MF2 = MF)
  
  chr_columns <- relations %>%
    map_lgl(is.character)
  
  relations[,chr_columns] <- relations[,chr_columns] %>%
    {
      .[is.na(.)] <- ''
      return(.)
    }
  
  relations <- relations %>%
    filter(Isotope2 == Isotope & Adduct2 == Adduct) %>%
    select(-(Isotope:`Measured m/z`))
  
  if (nrow(relations) > 0) {
    relations[relations == ''] <- NA  
  }
  
  return(relations)
}