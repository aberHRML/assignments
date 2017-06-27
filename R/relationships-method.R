#' @importFrom parallel makeCluster parApply stopCluster
#' @importFrom dplyr mutate bind_rows tbl_df filter
#' @importFrom dplyr inner_join semi_join select
#' @importFrom stringr str_sub str_replace_all
#' @importFrom mzAnnotation relationshipPredictor
#' @importFrom magrittr %>%
#' @importFrom tibble tibble

setMethod('relationships',signature = 'Annotation',
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
            clus <- makeCluster(parameters@nCores)
            rel <- parApply(clus,cors[1:100,c('Feature1','Feature2')],1,function(y,limit){
              mzAnnotation::relationshipPredictor(y,mode = c('n','p'),limit = limit)
            },limit = parameters@limit)
            stopCluster(clus)
            
            rel <- suppressWarnings(bind_rows(rel)) %>% 
              tbl_df() %>% 
              mutate(Error = as.numeric(Error)) %>%
              filter(Adduct1 %in% adducts$Adduct, Adduct2 %in% adducts$Adduct) %>%
              inner_join(cors,by = c('mz1' = 'Feature1','mz2' = 'Feature2')) %>%
              semi_join(adducts,by = c('Mode1' = 'Mode', 'Adduct1' = 'Adduct')) %>%
              semi_join(adducts,by = c('Mode2' = 'Mode', 'Adduct2' = 'Adduct')) %>% 
              filter(Isotope1 %in% c(NA,parameters@isotopes) & Isotope2 %in% c(NA,parameters@isotopes)) %>%
              select(mz1:r)
            
            x@relationships <- rel
            
            return(x)
          })