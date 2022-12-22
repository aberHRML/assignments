
elapsedTime <- function(start_time,end_time){
  {end_time - start_time} %>%
    .[3] %>%
    round(1) %>%
    seconds_to_period() %>%
    str_c('[',.,']')
}

#' @importFrom dplyr bind_cols all_of

eliminate <- function(MFs,by,direction){
  direct <- get(direction)
  
  MFs %>%
    bind_cols(MFs %>% select(by = by)) %>%
    group_by(Feature) %>% 
    filter(by == direct(by)) %>% 
    select(-all_of(by)) %>% 
    ungroup()
}

#' @importFrom dplyr rename
#' @importFrom purrr map_lgl

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
    group_split() %>% 
    furrr::future_map_dfr(~.x %>% 
                            mutate(M = calcM(mz,
                                             adduct = Adduct,
                                             isotope = Isotope))
                            ) %>% 
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

AIS <- function(assignment){
  possible_products <- expand_grid(
    Adduct = adducts(assignment) %>% 
      flatten_chr(),
    Isotope = c(NA,isotopes(assignment))
  )
  
  adducts_scores <- assignment %>% 
    adducts() %>% 
    map(~tibble(Adduct = .x,
                Adduct_Score = (length(.x) -1):0)) %>% 
    bind_rows()
  
  isotopes_scores <- assignment %>% 
    isotopes() %>% 
    {c(NA,.)} %>% 
    {tibble(Isotope = .,
            Isotope_Score = (length(.) -1):0)}
  
  possible_products %>% 
    left_join(adducts_scores,by = 'Adduct') %>% 
    left_join(isotopes_scores,by = 'Isotope') %>% 
    mutate(AIS = Adduct_Score + Isotope_Score,
           AIS = AIS / max(AIS)) %>% 
    select(-contains('Score'))
}

#' @importFrom tidyr expand_grid
#' @importFrom purrr flatten_chr map_dfr

maxAIS <- function(assignment){
  assignment_adducts <- adducts(assignment)
  assignment_isotopes <- isotopes(assignment)
  
  n_adducts <- assignment %>% 
    adducts() %>% 
    flatten_chr() %>% 
    length()
  
  n_isotopes <- assignment %>% 
    isotopes() %>% 
    {c(NA,.)} %>% 
    length()
  
  max_score <- (n_adducts * n_isotopes) / 2
  
  return(max_score)
}


#' @importFrom furrr future_map_dfr

generateMFs <- function(M,
                        ppm,
                        rank_threshold,
                        adduct_rules,
                        isotope_rules,
                        AIS){
  nM <- nrow(M)
  
  M %>%
    ungroup() %>%
    slice_sample(n = nM) %>%
    split(1:nrow(.)) %>%
    future_map_dfr(~{
      mf <- ipMF(mz = .x$mz,
                 adduct = .x$Adduct,
                 isotope = .x$Isotope,
                 ppm = ppm,
                 adduct_rules_table = adduct_rules,
                 isotope_rules_table = isotope_rules) %>% 
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
          left_join(AIS,
                    by = c('Adduct','Isotope'))
      } else {
        return(NULL)
      }
    },
    .options = furrr_options(seed = 1234))
}

#' @importFrom mzAnnotation transformationPossible

sanitiseTransformations <- function(graph_edges,transformation_rules_table){
  transforms <- dplyr::bind_rows(
    graph_edges %>% 
      filter(is.na(Transformation1)) %>% 
      select(
        from = MF1,
        to = MF2,
        transformation = Transformation2),
    graph_edges %>% 
      filter(is.na(Transformation2)) %>% 
      select(
        from = MF2,
        to = MF1,
        transformation = Transformation1)
    
  )
  
  transformation_possible <- transforms %>% 
    rowwise() %>% 
    group_split() %>% 
    map_lgl(~mzAnnotation::transformationPossible(
      .x$from,
      .x$to,
      .x$transformation,
      transformation_rules_table))
  
  graph_edges %>% 
    filter(transformation_possible)
}
