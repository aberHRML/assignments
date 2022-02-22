
setGeneric("calcRelationships", function(assignment,transformations = TRUE)
  standardGeneric("calcRelationships"))

#' @importFrom furrr future_map
#' @importFrom dplyr mutate bind_rows filter vars contains
#' @importFrom dplyr inner_join semi_join select mutate_at relocate
#' @importFrom stringr str_sub str_replace_all str_remove
#' @importFrom mzAnnotation relationshipCalculator
#' @importFrom magrittr %>%
#' @importFrom tibble tibble enframe

setMethod('calcRelationships',signature = 'Assignment',
          function(assignment,transformations = T){
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              message(blue('Calculating relationships '),cli::symbol$continue,'\r',appendLF = 'FALSE')
            }
            
            parameters <- as(assignment,'AssignmentParameters')
            
            cors <- assignment@preparedCorrelations
            
            if (isTRUE(transformations)) {
              trans <- c(NA,transformations(assignment))
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
                    setNames(str_remove(names(.),'1')),
                  .x %>% 
                    select(contains('2')) %>% 
                    setNames(str_remove(names(.),'2'))
                )
                
                modes <- mzs$Mode %>% 
                  unique()
                
                if (length(modes) > 1){
                  specified_adducts <- adducts(assignment) %>% 
                    unlist()
                } else {
                  specified_adducts <- adducts(assignment)[[modes]]
                }
                
                relationships <- relationshipCalculator(mzs$`m/z`,
                                                        limit = limit(assignment),
                                                        adducts = specified_adducts, 
                                                        isotopes = c(NA,isotopes(assignment)),
                                                        transformations = trans,
                                                        adduct_rules_table = adductRules(assignment),
                                                        isotope_rules_table = isotopeRules(assignment),
                                                        transformation_rules_table = transformationRules(assignment)) %>% 
                  left_join(mzs,by = c('m/z1' = 'm/z')) %>% 
                  rename(Mode1 = Mode) %>% 
                  left_join(mzs,by = c('m/z2' = 'm/z')) %>% 
                  rename(Mode2 = Mode) %>% 
                  relocate(contains('Mode'),.after = `m/z2`)
                
                if (length(modes) > 1){
                  adduct_modes <- adducts(assignment) %>% 
                    map(enframe,value = 'Adduct') %>% 
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
            
            relationships(assignment) <- rel
            
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