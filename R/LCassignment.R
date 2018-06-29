
LCassignment <- function(element){
  methods <- list(
    prepCorrelations = function(assignment){
      parameters <- assignment@parameters
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
      assignment %>%
        prepCorrelations()
    },
    relationships = function(assignment){
      assignment %>% 
        relationships()
    },
    adductIsotopeAssignment = function(assignment){
      assignment %>%
        addIsoAssign()
    }
  )
  
  if (!is.null(elements)) {
    method <- methods[[element]]
  } 
  return(method)
}