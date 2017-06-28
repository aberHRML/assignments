
addIsoScore <- function(add,iso,addRank,isoRank){
  add <- tibble(Adduct = add)
  iso <- tibble(Isotope = iso)
  addRank <- tibble(Adduct = unlist(addRank), 
                    Rank = unlist(lapply(addRank,function(x){
                      return(1:length(x))})))
  isoRank <- tibble(Isotope = c(NA,isoRank), Rank = 1:(length(isoRank) + 1))
  
  add <- inner_join(add, addRank,by = c('Adduct' = 'Adduct'))
  add <- sum(add$Rank)
  iso <- inner_join(iso, isoRank,by = c('Isotope' = 'Isotope'))
  iso <- sum(iso$Rank)
  return(add + iso)
}