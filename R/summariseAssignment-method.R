#' summariseAssignment-Assignment
#' @rdname summariseAssignment
#' @description Summarise features assigned to moleuclar formulas.
#' @param assignment S4 object of class Assignment
#' @importFrom dplyr desc
#' @export

setMethod('summariseAssignment',signature = 'Assignment',
          function(assignment){
            assigned <- assignment %>%
              assignments() %>%
              split(.$MF) %>%
              map(~{
                d <- .
                d$Isotope[is.na(d$Isotope)] <- ''
                d <- d %>%
                  mutate(IIP = str_c(Isotope,Adduct,sep = ' ')) %>%
                  arrange(`Measured m/z`)
                tibble(MF = d$MF[1],Features = str_c(d$Feature,collapse = '; '),`Isotopes & Ionisation Products` = str_c(d$IIP,collapse = '; '),Count = nrow(d))
              }) %>%
              bind_rows() %>%
              arrange(desc(Count))
            return(assigned)
})
