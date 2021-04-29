#' @importFrom furrr future_map
#' @importFrom dplyr mutate bind_rows filter vars
#' @importFrom dplyr inner_join semi_join select mutate_at
#' @importFrom stringr str_sub str_replace_all
#' @importFrom mzAnnotation relationshipCalculator
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom mzAnnotation adducts isotopes

setMethod('relationships',signature = 'Assignment',
          function(assignment,transformations = T){
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              message(blue('Calculating relationships '),cli::symbol$continue,'\r',appendLF = 'FALSE')
            }
            
            parameters <- assignment@parameters
            
            cors <- assignment@preparedCorrelations
            
            if (isTRUE(transformations)) {
              trans <- parameters@transformations
            } else {
              trans <- c()
            }
            
            rel <- cors %>%
              select(`m/z1`,`m/z2`,Mode1,Mode2) %>%
              split(1:nrow(.)) %>%
              future_map(~{
                relationshipCalculator(.x %>% 
                                         select(`m/z1`,`m/z2`) %>% 
                                         unlist(),
                                       limit = parameters@limit,
                                       modes = .x %>% 
                                         select(Mode1,Mode2) %>% 
                                         unlist(),
                                       adducts = parameters@adducts, 
                                       isotopes = parameters@isotopes,
                                       transformations = trans,
                                       adductTable = parameters@adductRules,
                                       isotopeTable = parameters@isotopeRules,
                                       transformationTable = parameters@transformationRules)
              }) %>%
              bind_rows() %>%
              inner_join(cors,by = c('m/z1' = 'm/z1','m/z2' = 'm/z2')) %>%
              select(Feature1:Mode2,
                     `m/z1`,
                     `m/z2`,
                     RetentionTime1,
                     RetentionTime2,
                     Adduct1:Transformation2,
                     log2IntensityRatio,
                     r,
                     Error,
                     ID) %>%
              mutate_at(vars(RetentionTime1,RetentionTime2),as.numeric)
            
            assignment@relationships <- rel
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              message(blue('Calculating relationships '),'\t',green(cli::symbol$tick),' ',elapsed)
            }
            
            return(assignment)
          })