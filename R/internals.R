#' @importFrom dplyr bind_cols

eliminate <- function(MFs,by,direction){
  direct <- get(direction)
  
  MFs %>%
    bind_cols(MFs %>% select(by = by)) %>%
    group_by(Feature) %>% 
    filter(by == direct(by)) %>% 
    select(-by) %>% 
    ungroup()
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

collateM <- function(rel,max_M){
  bind_rows(select(rel,
                   mz = `m/z1`,
                   RetentionTime = RetentionTime1,
                   Isotope = Isotope1,
                   Adduct = Adduct1, 
                   Feature = Feature1),
            select(rel,
                   mz = `m/z2`,
                   RetentionTime = RetentionTime2,
                   Isotope = Isotope2, 
                   Adduct = Adduct2, 
                   Feature = Feature2)) %>%
    distinct() %>%
    arrange(mz) %>%
    rowwise() %>%
    mutate(M = calcM(mz,
                     adduct = Adduct,
                     isotope = Isotope)) %>% 
    arrange(M) %>%
    filter(M <= max_M)
}

collateMFs <- function(rel,MF){
  bind_rows(select(rel,
                   Name = Name1,
                   Feature = Feature1,
                   mz = `m/z1`,
                   RetentionTime = RetentionTime1,
                   Isotope = Isotope1, 
                   Adduct = Adduct1, 
                   MF = MF1),
            select(rel,
                   Name = Name2,
                   Feature = Feature2,
                   mz = `m/z2`,
                   RetentionTime = RetentionTime2,
                   Isotope = Isotope2, 
                   Adduct = Adduct2,
                   MF = MF2)) %>%
    mutate(RetentionTime = as.numeric(RetentionTime)) %>%
    arrange(mz) %>%
    select(-mz) %>%
    left_join(MF, by = c("Feature", 
                         "RetentionTime",
                         "Isotope", 
                         "Adduct",
                         'MF')) %>%
    distinct() %>%
    mutate(ID = 1:nrow(.))
}

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

#' @importFrom tidyr expand_grid
#' @importFrom purrr flatten_chr

maxAddIsoScore <- function(assignment){
  
  assignment_adducts <- adducts(assignment)
  assignment_isotopes <- isotopes(assignment)
  
  adduct_scores <- assignment_adducts %>% 
    map(~tibble(adduct = .) %>% 
          mutate(adduct_rank = (nrow(.) - 1):0)) %>% 
    bind_rows(.id = 'mode')
  
  isotope_scores <- assignment_isotopes %>% 
    c('NA',.) %>% 
    tibble(isotope = .) %>% 
    mutate(isotope_rank = (nrow(.) - 1):0)
  
  max_total <- assignment_adducts %>% 
    map_dfr(~tibble(max_total = length(.x) - 1),
            .id = 'mode') %>% 
    mutate(max_total = max_total + 
             length(assignment_isotopes))
  
  add_iso_combs <- expand_grid(adduct = assignment_adducts %>% 
                                 flatten_chr(),
                               isotope = c('NA',assignment_isotopes)
  )
  
  add_iso_scores <- add_iso_combs %>% 
    left_join(adduct_scores, 
              by = "adduct") %>% 
    left_join(isotope_scores, 
              by = "isotope") %>% 
    mutate(total = adduct_rank + 
             isotope_rank) %>% 
    left_join(max_total,by = 'mode') %>% 
    mutate(score = total / max_total)
  
  max_score <- sum(add_iso_scores$score)
  
  return(max_score)
}


#' @importFrom furrr future_map_dfr

generateMFs <- function(M,
                        ppm,
                        rank_threshold,
                        adducts,
                        isotopes){
  nM <- nrow(M)
  
  M %>%
    ungroup() %>%
    slice_sample(n = nM) %>%
    split(1:nrow(.)) %>%
    future_map_dfr(~{
      mf <- ipMF(mz = .x$mz,
                 adduct = .x$Adduct,
                 isotope = .x$Isotope,
                 ppm = ppm) %>% 
        mutate(Rank = rank(100 - `Plausibility (%)`,
                           ties.method = 'min')) %>% 
        filter(Rank <= rank_threshold)
      
      if (nrow(mf) > 0) {
        mf %>%
          left_join(select(M,
                           Feature,
                           RetentionTime,
                           M,
                           mz),
                    by = c('Measured M' = 'M','Measured m/z' = 'mz')) %>% 
          rowwise() %>%
          select(Feature,RetentionTime,MF,Isotope,Adduct,`Theoretical M`,
                 `Measured M`,`Theoretical m/z`,`Measured m/z`, `PPM error`,
                 `MF Plausibility (%)` = `Plausibility (%)`) %>%
          
          rowwise() %>%
          mutate(AddIsoScore = addIsoScore(Adduct,
                                           Isotope,
                                           adducts(assignment),
                                           isotopes(assignment))) %>%
          ungroup()
      } else {
        return(NULL)
      }
    },
    .options = furrr_options(seed = 1234))
}
