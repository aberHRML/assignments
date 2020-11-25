
setMethod('prepCorrelations',signature = 'Assignment',
          function(assignment){
            
            if (assignment@log$verbose == T) {
             startTime <- proc.time()
             message(blue('Preparing correlations '),cli::symbol$continue,'\r',appendLF = FALSE)
            }
            
            correlations <- assignment@preparedCorrelations
            
            correlations <- correlations %>%
              mutate(Mode1 = str_split_fixed(Feature1,'@',2) %>% 
                       .[,1] %>%
                       str_sub(1,1),
                     Mode2 = str_split_fixed(Feature2,'@',2) %>% 
                       .[,1] %>%
                       str_sub(1,1),
                     `m/z1` = str_split_fixed(Feature1,'@',2) %>% 
                       .[,1] %>% 
                       str_replace_all('[:alpha:]','') %>% 
                       as.numeric(),
                     `m/z2` = str_split_fixed(Feature2,'@',2) %>% 
                       .[,1] %>% 
                       str_replace_all('[:alpha:]','') %>% 
                       as.numeric(),
                     RetentionTime1 = str_split_fixed(Feature1,'@',2) %>% 
                       .[,2] %>%
                       as.numeric(),
                     RetentionTime2 = str_split_fixed(Feature2,'@',2) %>% 
                       .[,2] %>%
                       as.numeric(),
                     ID = 1:nrow(.)
              ) %>%
              select(Feature1,Feature2,Mode1:RetentionTime2,log2IntensityRatio,r,ID)
            
            assignment@preparedCorrelations <- correlations
            
            if (assignment@log$verbose == T) {
              endTime <- proc.time()
              elapsed <- {endTime - startTime} %>%
                .[3] %>%
                round(1) %>%
                seconds_to_period() %>%
                str_c('[',.,']')
              message(blue('Preparing correlations '),'\t\t',green(cli::symbol$tick),' ',elapsed)
            }
            
            return(assignment)
          })