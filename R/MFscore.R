#' @importFrom CHNOSZ count.elements

MFscore <- function(mf){
  elements <- c('C','H','N','O','P','S')
  eleFreq <- as.vector(count.elements(mf))
  ele <- names(count.elements(mf))
  
  if (length(which(!(elements %in% ele))) > 0) {
    eleFreq <- c(eleFreq,rep(0,length(which(!(elements %in% ele)))))
    ele <- c(ele,elements[!(elements %in% ele)])
  }
  
  colnames(eleFreq) <- NULL
  
  eleRatios <- c(
    `H/C` = if ('H' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'H']/eleFreq[ele == 'C']
    },
    `N/C` = if ('N' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'N']/eleFreq[ele == 'C']
    },
    `O/C` = if ('O' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'O']/eleFreq[ele == 'C']
    },
    `P/C` = if ('P' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'P']/eleFreq[ele == 'C']
    },
    `S/C` = if ('S' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'S']/eleFreq[ele == 'C']
    },
    `N/O` = if ('N' %in% ele & 'O' %in% ele) {
      eleFreq[ele == 'N']/eleFreq[ele == 'O']
    },
    `P/O` = if ('P' %in% ele & 'O' %in% ele) {
      eleFreq[ele == 'P']/eleFreq[ele == 'O']
    },
    `S/O` = if ('S' %in% ele & 'O' %in% ele) {
      eleFreq[ele == 'S']/eleFreq[ele == 'O']
    },
    `O/P` = if ('O' %in% ele & 'P' %in% ele) {
      eleFreq[ele == 'O']/eleFreq[ele == 'P']
    },
    `S/P` = if ('S' %in% ele & 'P' %in% ele) {
      eleFreq[ele == 'S']/eleFreq[ele == 'P']
    }
  )
  if (!is.null(eleRatios)) {
    if ('H/C' %in% names(eleRatios)) {
      eleRatios['H/C'] <- abs(eleRatios['H/C'] - 1.6)
    }
    if ('O/C' %in% names(eleRatios)) {
      eleRatios['O/C'] <- abs(eleRatios['O/C'] - 0.3)
    }
    if  (is.nan(eleRatios['N/O'])) {
      eleRatios['N/O'] <- 0
    }
    if (is.infinite(eleRatios['N/O'])) {
      eleRatios['N/O'] <- eleFreq[ele == 'N']
    }
    if (is.nan(eleRatios['P/O'])) {
      eleRatios['P/O'] <- 0
    }
    if (is.infinite(eleRatios['P/O'])) {
      eleRatios['P/O'] <- eleFreq[ele == 'P']
    }
    if (is.nan(eleRatios['S/O'])) {
      eleRatios['S/O'] <- 0
    }
    if (is.infinite(eleRatios['S/O'])) {
      eleRatios['S/O'] <- eleFreq[ele == 'S']
    }
    if (is.nan(eleRatios['O/P'])) {
      eleRatios['O/P'] <- 0
    }
    if ('O/P' %in% names(eleRatios) & eleRatios['O/P'] >= 3) {
      eleRatios['O/P'] <- 0
    }
    if (is.nan(eleRatios['P/S'])) {
      eleRatios['P/S'] <- 0
    }
    if (is.nan(eleRatios['S/P'])) {
      eleRatios['S/P'] <- 0
    }
    if (is.infinite(eleRatios['S/P'])) {
      eleRatios['S/P'] <- eleFreq[ele == 'S']
    }
    score <- sum(eleRatios) 
  } else {
    score <- NA
  }
  
  return(score)
}
