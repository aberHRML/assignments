
setMethod('summariseAssignment',signature = 'Assignment',
          function(assignment){
            assigned <- assignment %>%
              assignments() %>%
              split(.$MF) %>%
              map(~{
                d <- .
                d$Isotope[is.na(d$Isotope)] <- ''
                d <- d %>%
                  mutate(Feature = str_c(Mode,`Measured m/z`),IIP = str_c(Isotope,Adduct,sep = ' ')) %>%
                  arrange(`Measured m/z`)
                tibble(MF = d$MF[1],Features = str_c(d$Feature,collapse = '; '),IIPs = str_c(d$IIP,collapse = '; '),Count = nrow(d))
              }) %>%
              bind_rows() %>%
              arrange(desc(Count))
            return(assigned)
})
