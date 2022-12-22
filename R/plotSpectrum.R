
spectrumPlot <- function(dat,MF){
  check_installed(c('ggplot2',
                    'ggrepel',
                    'ggthemes'))
  
  dat$Mode[dat$Mode == 'n'] <- 'Negative mode'
  dat$Mode[dat$Mode == 'n'] <- 'Positive mode'
  
  ggplot2::ggplot(dat) +
    ggplot2::geom_segment(
      ggplot2::aes(x = `m/z`,
                   xend = `m/z`, 
                   y = 0, 
                   yend = `Relative Abundance`),
      colour = ggthemes::ptol_pal()(1)) +
    ggrepel::geom_text_repel(
      ggplot2::aes(x = `m/z`,
                   y = `Relative Abundance`,
                   label = Label)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(),
          axis.title = ggplot2::element_text(face = 'bold'),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = 'bold'),
          plot.title = ggplot2::element_text(face = 'bold',
                                    hjust = 0.5)) +
    ggplot2::labs(title = MF,
         y = 'Relative Abundance') +
    ggplot2::facet_wrap(~Mode,scales = 'free')
}

#' plotSpectrum
#' @rdname plotSpectrum
#' @description Plot a spectrum for a given molecular formula 
#' @param assignment S4 object of class Assignment
#' @param MF molecular formula
#' @importFrom tidyr gather
#' @importFrom dplyr group_by summarise
#' @export

setGeneric('plotSpectrum',function(assignment,MF)
  standardGeneric('plotSpectrum'))

#' @rdname plotSpectrum

setMethod('plotSpectrum',signature = 'Assignment',function(assignment,MF){
  mf <- MF
  
  feat <- assignment %>%
    assignments() %>%
    filter(MF == mf)
  
  dat <- assignment@data %>%
    select(all_of(feat$Feature)) %>%
    gather('Feature','Intensity') %>%
    group_by(Feature) %>%
    summarise(Intensity = mean(Intensity)) %>%
    mutate(`Relative Abundance` = Intensity / max(Intensity)) %>%
    left_join(feat %>%
                select(Feature,Adduct,Isotope,Mode,`m/z` = `Measured m/z`), by = "Feature") 
  
  dat[is.na(dat)] <- ''
  
  dat <- dat %>%
    mutate(Label = str_c(Isotope,Adduct,sep = ' '))
  
  spectrumPlot(dat,MF)
})