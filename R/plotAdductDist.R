
plotDist <- function(x){
  check_installed(c('ggplot2',
                    'ggthemes'))
  
  ggplot2::ggplot(x,ggplot2::aes(x = Adduct)) + 
    ggplot2::geom_bar(colour = 'black',
                      fill = ggthemes::ptol_pal()(1)) + 
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::theme_bw() + 
    ggplot2::facet_wrap(~Isotope,
               scales = 'free') + 
    ggplot2::labs(title = x$Mode[1],
         y = 'Count',
         caption = str_c('N = ',nrow(x))) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = 'bold',hjust = 0.5),
          axis.title = ggplot2::element_text(face = 'bold'),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(),
          panel.grid = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = 'bold'))
}

#' @rdname plotting
#' @importFrom tidyr replace_na
#' @export

setGeneric('plotAdductDist',function(assignment){
  standardGeneric('plotAdductDist')
})

#' @rdname plotting
#' @importFrom rlang check_installed

setMethod('plotAdductDist',signature = 'Assignment',
          function(assignment){
            
            check_installed('patchwork')
            
            assign <- assignment %>%
              assignments() %>% 
              replace_na(list(Isotope = '')) %>% 
              mutate(Isotope = factor(Isotope,
                                      levels = c('',isotopes(assignment)))
                     )
            
            assign$Mode[assign$Mode == 'n'] <- 'Negative Mode'
            assign$Mode[assign$Mode == 'p'] <- 'Positive Mode'
            
            assign %>% 
              split(.$Mode) %>% 
              map(plotDist) %>% 
              patchwork::wrap_plots()
          }
)