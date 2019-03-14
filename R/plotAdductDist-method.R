#' plotAdductDist-Assignment
#' @importFrom patchwork wrap_plots
#' @importFrom ggthemes ptol_pal
#' @importFrom ggplot2 ggplot geom_bar theme_bw facet_wrap theme element_text
#' @export

setMethod('plotAdductDist',signature = 'Assignment',
          function(assignment){
            assign <- assignment %>%
              assignments()
            
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
                  theme(plot.title = element_text(face = 'bold'),
                        axis.title = element_text(face = 'bold'),
                        axis.text.x = element_text(angle = 45, hjust = 1))
                }) %>% 
              wrap_plots()
          }
)