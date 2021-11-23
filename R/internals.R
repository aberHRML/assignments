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

#' @importFrom dplyr bind_cols

eliminate <- function(MFs,by,direction){
  MFs %>%
    bind_cols(MFs %>% select(by = by)) %>%
    split(.$Feature) %>%
    map(~{
      d <- .
      direct <- get(direction)
      d %>%
        filter(by == direct(by))
    }) %>%
    bind_rows() %>%
    select(-by)
}

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
      .
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
      .
    }
  
  relations <- relations %>%
    filter(Isotope2 == Isotope & Adduct2 == Adduct) %>%
    select(-(Isotope:`Measured m/z`))
  
  if (nrow(relations) > 0) {
    relations[relations == ''] <- NA  
  }
  
  return(relations)
}

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