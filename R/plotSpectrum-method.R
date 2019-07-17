#' plotSpectrum
#' @rdname plotSpectrum
#' @description Plot a spectrum for a given molecular formula 
#' @param assignment S4 object of class Assignment
#' @param MF molecular formula
#' @importFrom tidyr gather
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 geom_segment
#' @importFrom ggrepel geom_text_repel
#' @export

setMethod('plotSpectrum',signature = 'Assignment',function(assignment,MF){
  mf <- MF
  
  feat <- assignment %>%
    assignments() %>%
    filter(MF == mf)
  
  dat <- assignment@data %>%
    select(feat$Feature) %>%
    gather('Feature','Intensity') %>%
    group_by(Feature) %>%
    summarise(Intensity = mean(Intensity)) %>%
    mutate(`Relative Abundance` = Intensity / max(Intensity)) %>%
    left_join(feat %>%
                select(Feature,Adduct,Isotope,Mode,`m/z` = `Measured m/z`), by = "Feature") 
  
  dat[is.na(dat)] <- ''
  
  dat <- dat %>%
    mutate(Label = str_c(Isotope,Adduct,sep = ' '))
  
  ggplot(dat) +
    geom_segment(aes(x = `m/z`,xend = `m/z`, y = 0, yend = `Relative Abundance`),colour = ptol_pal()(1)) +
    geom_text_repel(aes(x = `m/z`,y = `Relative Abundance`,label = Label)) +
    theme_bw() +
    labs(title = MF,
         y = 'Relative Abundance') +
    facet_wrap(~Mode,scales = 'free')
})