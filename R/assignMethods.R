
assignMethods <- function(method = NULL) {
  methods <- list(
    FIE = function(assignment) {
      assignment <- assignment %>% 
        prepCorrelations() %>%
        relationships() %>% 
        addIsoAssign()
      count <- 0
      while (T) {
        count <- count + 1
        assignment <- suppressWarnings(transformationAssign(assignment))
        if (nrow(assignment@transAssign[[count]]$assigned) == 0) {
          assignment@transAssign <- assignment@transAssign[-count] 
          break()
        }
      }
    `RP-LC` = function(assignment){
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
    }
      return(assignment)
    }
  )
  
  if (is.null(method)) {
    method <- methods
  } else {
    method <- methods[[method]]
  }
  return(method)
}