
setMethod('assignments',signature = 'Assignment',
          function(assignment){
            assignment@assignments %>%
              mutate(Feature = str_c(Mode,`Measured m/z`)) %>%
              select(Feature,RetentionTime:Mode)
})
