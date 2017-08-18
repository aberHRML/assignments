#' @importFrom parallel makeCluster parApply stopCluster
#' @importFrom dplyr mutate bind_rows tbl_df filter
#' @importFrom dplyr inner_join semi_join select
#' @importFrom stringr str_sub str_replace_all
#' @importFrom mzAnnotation relationshipPredictor
#' @importFrom magrittr %>%
#' @importFrom tibble tibble

setMethod('relationships',signature = 'Assignment',
          function(x){
            parameters <- x@parameters
            
            adducts <- lapply(parameters@adducts,function(y){tibble(Adduct = y)})
            adducts <- bind_rows(adducts,.id = 'Mode')
            
            cors <- x@correlations %>%
              mutate(Mode1 = str_sub(Feature1,1,1),
                     Mode2 = str_sub(Feature2,1,1),
                     Feature1 = as.numeric(str_replace_all(Feature1,'[:alpha:]','')), 
                     Feature2 = as.numeric(str_replace_all(Feature2,'[:alpha:]',''))
              )
            clus <- makeCluster(parameters@nCores,type = parameters@clusterType)
            rel <- parApply(clus,select(cors,Feature1,Feature2,Mode1,Mode2),1,function(y,limit,add,iso,trans){
              add <- unlist(add[names(add) %in% y[3:4]])
              mzAnnotation::relationshipPredictor(as.numeric(y[1:2]),limit = limit,adducts = add,isotopes = iso,transformations = trans)
            },limit = parameters@limit, add = parameters@adducts, iso = parameters@isotopes,trans = parameters@transformations)
            stopCluster(clus)
            
            rel <- suppressWarnings(bind_rows(rel)) %>% 
              tbl_df() %>% 
              mutate(Error = as.numeric(Error)) %>%
              filter(Adduct1 %in% adducts$Adduct, Adduct2 %in% adducts$Adduct) %>%
              inner_join(cors,by = c('m/z1' = 'Feature1','m/z2' = 'Feature2')) %>%
              semi_join(adducts,by = c('Mode1' = 'Mode', 'Adduct1' = 'Adduct')) %>%
              semi_join(adducts,by = c('Mode2' = 'Mode', 'Adduct2' = 'Adduct')) %>% 
              filter(Isotope1 %in% c(NA,parameters@isotopes) & Isotope2 %in% c(NA,parameters@isotopes)) %>%
              select(`m/z1`:r)
            
            x@relationships <- rel
            
            return(x)
          })