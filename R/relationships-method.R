#' @importFrom parallel makeCluster parApply stopCluster
#' @importFrom dplyr mutate bind_rows tbl_df filter
#' @importFrom dplyr inner_join semi_join select
#' @importFrom stringr str_sub str_replace_all
#' @importFrom mzAnnotation relationshipPredictor
#' @importFrom magrittr %>%
#' @importFrom tibble tibble

setMethod('relationships',signature = 'Assignment',
          function(assignment){
            parameters <- assignment@parameters
            
            cors <- assignment@correlations
            
            clus <- makeCluster(parameters@nCores,type = parameters@clusterType)
            rel <- parApply(clus,select(cors,`m/z1`,`m/z2`,Mode1,Mode2),1,function(y,limit,add,iso,trans){
              mzAnnotation::relationshipPredictor(as.numeric(y[1:2]),limit = limit,modes = y[3:4],adducts = add,isotopes = iso,transformations = trans)
            },limit = parameters@limit, add = parameters@adducts, iso = parameters@isotopes,trans = parameters@transformations) %>%
              bind_rows() %>%
              as_tibble() %>%
              inner_join(cors,by = c('m/z1' = 'm/z1','m/z2' = 'm/z2')) %>%
              select(Feature1:Mode2,`m/z1`,`m/z2`,rt1,rt2,Adduct1:Transformation2,log2IntensityRatio,r,Error)
              
            stopCluster(clus)
            
            assignment@relationships <- rel
            
            return(assignment)
          })