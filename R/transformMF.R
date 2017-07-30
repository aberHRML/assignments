#' @importFrom tidyr gather
#' @importFrom stringr str_c str_replace

transformMF <- function(MF,transformation){
  if (!is.na(transformation)) {
    elements <- c('C','H','O','N','P','S')
    
    transformations <- rules$transformation
    transformation <- filter(transformations,`MF Change` == transformation) %>%
      select(C:S)
    
    MF <- makeup(MF)
    
    if (length(which(!(elements %in% names(MF)))) > 0) {
      MF <- c(MF,rep(0,length(which(!(elements %in% names(MF))))))
      names(MF)[names(MF) == ''] <- elements[!(elements %in% names(MF))]
    }
    
    MF <- MF[order(names(MF))]
    MF <- MF + transformation
    MF <- gather(MF,'Element','Frequency') 
    
    if (T %in% (MF$Frequency < 0)) {
      MF <- NA
    } else {
      MF <- MF %>%
        filter(Frequency > 0)
      
      MF$Frequency[MF$Frequency == 1] <- ''
      
      MF <- str_c(MF$Element,MF$Frequency) %>%
        str_c(collapse = '') 
    }
  }
  return(MF)
}