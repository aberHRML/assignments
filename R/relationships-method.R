#' @importFrom furrr future_map
#' @importFrom dplyr mutate bind_rows filter vars contains
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
            
            parameters <- as(assignment,'AssignmentParameters')
            
            cors <- assignment@preparedCorrelations
            
            if (isTRUE(transformations)) {
              trans <- c(NA,parameters@transformations)
            } else {
              trans <- NA
            }
            
            rel <- cors %>%
              select(`m/z1`,`m/z2`,Mode1,Mode2) %>%
              split(1:nrow(.)) %>%
              future_map(~{
                
                mzs <- bind_rows(
                  .x %>% 
                    select(contains('1')) %>% 
                    setNames(stringr::str_remove(names(.),'1')),
                  .x %>% 
                    select(contains('2')) %>% 
                    setNames(stringr::str_remove(names(.),'2'))
                )
                
                modes <- mzs$Mode %>% 
                  unique()
                
                if (length(modes) > 1){
                  adducts <- parameters@adducts %>% 
                    unlist()
                } else {
                  adducts <- parameters@adducts[[modes]]
                }
                
                relationships <- relationshipCalculator(mzs$`m/z`,
                                                        limit = parameters@limit,
                                                        adducts = adducts, 
                                                        isotopes = c(NA,parameters@isotopes),
                                                        transformations = trans,
                                                        adductTable = parameters@adductRules,
                                                        isotopeTable = parameters@isotopeRules,
                                                        transformationTable = parameters@transformationRules) %>% 
                  left_join(mzs,by = c('m/z1' = 'm/z')) %>% 
                  rename(Mode1 = Mode) %>% 
                  left_join(mzs,by = c('m/z2' = 'm/z')) %>% 
                  rename(Mode2 = Mode) %>% 
                  dplyr::relocate(contains('Mode'),.after = `m/z2`)
                
                if (length(modes) > 1){
                  adduct_modes <- parameters@adducts %>% 
                    map(tibble::enframe,value = 'Adduct') %>% 
                    bind_rows(.id = 'Mode') %>% 
                    select(-name)
                  
                  relationships <- relationships %>% 
                    inner_join(adduct_modes,
                               by = c('Mode1' = 'Mode',
                                      'Adduct1' = 'Adduct')) %>% 
                    inner_join(adduct_modes,
                               by = c('Mode2' = 'Mode',
                                      'Adduct2' = 'Adduct'))
                }
                
                return(relationships)
              }) %>%
              bind_rows() %>%
              inner_join(cors,by = c('m/z1','m/z2','Mode1','Mode2')) %>%
              select(contains('Feature'),
                     contains('Mode'),
                     contains('m/z'),
                     contains('RetentionTime'),
                     contains('Adduct'),
                     contains('Isotope'),
                     contains('Transformation'),
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