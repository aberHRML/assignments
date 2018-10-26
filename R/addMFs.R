#' @importFrom dplyr rename

addMFs <- function(rel,MF){

  relations <- rel %>%
    filter(Feature1 %in% MF$Feature,Feature2 %in% MF$Feature) %>%
    left_join(MF %>%
                select(Feature,MF,Isotope,Adduct,`Measured m/z`),by = c('Feature1' = 'Feature')) %>%
    rename(MF1 = MF)
  relations[is.na(relations)] <- ''
  
  relations <- relations %>%
    filter(Isotope1 == Isotope & Adduct1 == Adduct) %>%
    select(-(Isotope:`Measured m/z`)) %>%
    left_join(MF %>%
                select(Feature,MF,Isotope,Adduct,`Measured m/z`),by = c('Feature2' = 'Feature')) %>%
    rename(MF2 = MF)
  relations[is.na(relations)] <- ''
  
  relations <- relations %>%
    filter(Isotope2 == Isotope & Adduct2 == Adduct) %>%
    select(-(Isotope:`Measured m/z`))
  relations[relations == ''] <- NA
  
  return(relations)
}