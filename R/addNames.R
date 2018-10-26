
addNames <- function(rel){
  iso <- rel
  iso$Isotope1[is.na(iso$Isotope1)] <- ''
  iso$Isotope2[is.na(iso$Isotope2)] <- ''
  iso <- iso %>%
    mutate(Name1 = str_c(Feature1,MF1,Isotope1,Adduct1,sep = ' '),
           Name2 = str_c(Feature2,MF2,Isotope2,Adduct2,sep = ' '))
  rel %>%
    bind_cols(iso %>%
                select(Name1,Name2)) %>%
    select(Name1,Name2,Feature1:MF2)
}