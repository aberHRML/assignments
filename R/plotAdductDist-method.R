#' Plot adduct frequency distributions
#' @rdname plotAdductDist
#' @description Plot adduct frequency distributions.
#' @param assignment S4 object of class Assignment
#' @importFrom patchwork wrap_plots
#' @importFrom ggthemes ptol_pal
#' @importFrom ggplot2 ggplot geom_bar theme_bw facet_wrap theme element_text
#' @importFrom tidyr replace_na
#' @export

setGeneric('plotAdductDist',function(assignment){
  standardGeneric('plotAdductDist')
})

#' @rdname plotAdductDist

setMethod('plotAdductDist',signature = 'Assignment',
          function(assignment){
            
            isotopes <- assignmen
            
            assign <- assignment %>%
              assignments() %>% 
              replace_na(list(Isotope = '')) %>% 
              mutate(Isotope = factor(Isotope,levels = c()))
            
            assign$Mode[assign$Mode == 'n'] <- 'Negative Mode'
            assign$Mode[assign$Mode == 'p'] <- 'Positive Mode'
            
            assign %>% split(.$Mode) %>% 
              map(~{
                d <- .
                ggplot(d,aes(x = Adduct)) + 
                  geom_bar(colour = 'black',fill = ptol_pal()(1)) + 
                  theme_bw() + 
                  facet_wrap(~Isotope) + 
                  labs(title = d$Mode[1],
                       y = 'Count',
                       caption = str_c('N = ',nrow(d))) +
                  theme(plot.title = element_text(face = 'bold',hjust = 0.5),
                        axis.title = element_text(face = 'bold'),
                        axis.text.x = element_text(angle = 45, hjust = 1),
                        panel.border = element_blank(),
                        axis.line = element_line(),
                        panel.grid = element_blank(),
                        strip.background = element_blank(),
                        strip.text = element_text(face = 'bold'))
                }) %>% 
              wrap_plots()
          }
)