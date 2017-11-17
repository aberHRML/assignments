
LCassignment <- function(assignment){
    assignment@correlations <- assignment@correlations %>%
      mutate( rt1 = str_split_fixed(Feature1,'@',2) %>% 
                as_tibble() %>% 
                select(V2) %>%
                unlist() %>%
                as.numeric(),
              rt2 = str_split_fixed(Feature2,'@',2) %>% 
                as_tibble() %>% 
                select(V2) %>%
                unlist() %>%
                as.numeric(),
              rtDiff = abs(rt1 - rt2)
      ) %>%
      filter(rtDiff <= parameters@RTwindow) %>%
      select(Feature1:r)
    assignment <- assignment %>%
      prepCorrelations() %>%
      relationships() %>%
      addIsoAssign()
    
  return(assignment)
}