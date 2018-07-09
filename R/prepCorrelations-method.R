
setMethod('prepCorrelations',signature = 'Assignment',
          function(assignment){
            
            if (assignment@log$verbose == T) {
             startTime <- proc.time()
             cat(blue('Preparing correlations '),cli::symbol$continue,'\r',sep = '')
            }
            
            correlations <- assignment@correlations
            
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
                     RetentionTime1 = str_split_fixed(Feature1,'@',2) %>% 
                       as_tibble() %>% 
                       select(V2) %>%
                       unlist() %>%
                       as.numeric(),
                     RetentionTime2 = str_split_fixed(Feature2,'@',2) %>% 
                       as_tibble() %>% 
                       select(V2) %>%
                       unlist() %>%
                       as.numeric()
              ) %>%
              select(Feature1,Feature2,Mode1:RetentionTime2,log2IntensityRatio,r)
            
            assignment@correlations <- correlations
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              cat(blue('Preparing correlations '),'\t\t',green(cli::symbol$tick),' ',elapsed,'\n',sep = '')
            }
            
            return(assignment)
          })