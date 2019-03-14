#' @importFrom purrr map_lgl

addIsoScore <- function(add,iso,addRank,isoRank){
  add <- tibble(Adduct = add)
  iso <- tibble(Isotope = iso)
  iso$Isotope[is.na(iso$Isotope)] <- 'NA'
  addRank <- addRank[map_lgl(addRank,~{add %in% .})] %>%
    .[[1]] %>%
    {tibble(Adduct = ., 
                    Rank = (length(.) - 1):0)}
  isoRank <- tibble(Isotope = c('NA',isoRank), Rank = length(isoRank):0)
  
  add <- left_join(add, addRank,by = 'Adduct') %>%
    .$Rank
  iso <- left_join(iso, isoRank,by = 'Isotope') %>%
    .$Rank
  
  maxScore <- max(addRank$Rank) + max(isoRank$Rank)
  score <- (add + iso)/maxScore
  
  return(score)
}
