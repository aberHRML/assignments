#' assignMFs
#' @description assign molecular formulas to a set of given m/z.
#' @param correlations table containing correlations of m/z to assign molecular formulas
#' @param parameters an S4 object of class AssignmentParamters containing the parameters for molecular formula assignment
#' @importFrom tibble tibble
#' @importFrom stringr str_split_fixed
#' @examples 
#' \dontrun{
#' res <- assignMFs(correlations,assignmentParameters('FIE'))
#' }
#' @export

assignMFs <- function(correlations,parameters) {
  options(digits = 10)
  
  correlations <- correlations %>%
    mutate(Mode1 = str_split_fixed(Feature1,'@',2) %>% 
             as_tibble() %>% 
             select(V1) %>% 
             unlist() %>% 
             str_sub(1,1),
           Mode2 = str_split_fixed(Feature2,'@',2) %>% 
             as_tibble() %>% 
             select(V1) %>% 
             unlist() %>% 
             str_sub(1,1),
           `m/z1` = str_split_fixed(Feature1,'@',2) %>% 
             as_tibble() %>% 
             select(V1) %>% 
             unlist() %>% 
             str_replace_all('[:alpha:]','') %>% 
             as.numeric(),
           `m/z2` = str_split_fixed(Feature2,'@',2) %>% 
             as_tibble() %>% 
             select(V1) %>% 
             unlist() %>% 
             str_replace_all('[:alpha:]','') %>% 
             as.numeric(),
           rt1 = str_split_fixed(Feature1,'@',2) %>% 
             as_tibble() %>% 
             select(V2) %>%
             unlist() %>%
             as.numeric(),
           rt2 = str_split_fixed(Feature1,'@',2) %>% 
             as_tibble() %>% 
             select(V2) %>%
             unlist() %>%
             as.numeric()
    ) %>%
    select(Feature1,Feature2,Mode1:rt2,log2IntensityRatio,r)
  
  assignment <- new('Assignment',
                    parameters = parameters,
                    correlations = correlations,
                    relationships = tibble(),
                    addIsoAssign = list(),
                    transAssign = list(),
                    assignments  = tibble()
  )
  
  method <- assignMethods(assignment@parameters@technique)
  assignment <- method(assignment)
  return(assignment)
}