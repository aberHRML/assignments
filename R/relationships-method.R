#' @importFrom parallel makeCluster parApply stopCluster
#' @importFrom dplyr mutate bind_rows tbl_df filter
#' @importFrom dplyr inner_join semi_join select
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
            
            clus <- makeCluster(parameters@nCores,type = parameters@clusterType)
            rel <- cors %>%
              select(`m/z1`,`m/z2`,Mode1,Mode2) %>%
              split(1:nrow(.)) %>%
              parLapply(cl = clus,function(y,limit,add,iso,trans,addRules,isoRules,transRules){
              mzAnnotation::relationshipCalculator(y %>%
                                                     select(`m/z1`,`m/z2`) %>%
                                                     unlist(),
                                                   limit = limit,
                                                   modes = y %>%
                                                     select(Mode1,Mode2) %>%
                                                     unlist(),
                                                   adducts = add,
                                                   isotopes = iso,
                                                   transformations = trans,
                                                   adductTable = addRules,
                                                   isotopeTable = isoRules,
                                                   transformationTable = transRules)
            },
            limit = parameters@limit,
            add = parameters@adducts, 
            iso = parameters@isotopes,
            trans = trans,
            addRules = parameters@adductRules,
            isoRules = parameters@isotopeRules,
            transRules = parameters@transformationRules
            ) %>%
              bind_rows() %>%
              inner_join(cors,by = c('m/z1' = 'm/z1','m/z2' = 'm/z2')) %>%
              select(Feature1:Mode2,`m/z1`,`m/z2`,RetentionTime1,RetentionTime2,Adduct1:Transformation2,log2IntensityRatio,r,Error,ID) %>%
              mutate(RetentionTime1 = as.numeric(RetentionTime1),RetentionTime2 = as.numeric(RetentionTime2))
            
            stopCluster(clus)
            
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