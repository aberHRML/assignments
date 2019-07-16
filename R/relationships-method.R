#' @importFrom parallel makeCluster parApply stopCluster
#' @importFrom dplyr mutate bind_rows tbl_df filter
#' @importFrom dplyr inner_join semi_join select
#' @importFrom stringr str_sub str_replace_all
#' @importFrom mzAnnotation relationshipCalculator
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom mzAnnotation adducts isotopes

setMethod('relationships',signature = 'Assignment',
          function(assignment){
            
            if (assignment@log$verbose == T) {
              startTime <- proc.time()
              cat(blue('Calculating relationships '),cli::symbol$continue,'\r',sep = '')
            }
            
            parameters <- assignment@parameters
            
            cors <- assignment@correlations
            
            clus <- makeCluster(parameters@nCores,type = parameters@clusterType)
            rel <- parApply(clus,select(cors,`m/z1`,`m/z2`,Mode1,Mode2),1,function(y,limit,add,iso,trans){
              mzAnnotation::relationshipCalculator(as.numeric(y[1:2]),
                                                  limit = limit,
                                                  modes = y[3:4],
                                                  adducts = add,
                                                  isotopes = iso,
                                                  transformations = trans,
                                                  adductTable = adducts(),
                                                  isotopeTable = isotopes(),
                                                  transformationTable = transformations())
            },limit = parameters@limit, add = parameters@adducts, iso = parameters@isotopes,trans = parameters@transformations) %>%
              bind_rows() %>%
              as_tibble() %>%
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
              cat(blue('Calculating relationships '),'\t',green(cli::symbol$tick),' ',elapsed,'\n',sep = '')
            }
            
            return(assignment)
          })