#' @importFrom CHNOSZ makeup

MFscore <- function(mf){
  eleFreq <- makeup(mf)

  eleRatios <- list(
    `H/C` = if ('H' %in% names(eleFreq) & 'C' %in% names(eleFreq)) {
      round(eleFreq['H']/eleFreq['C'],2)
    },
    `N/C` = if ('N' %in% names(eleFreq) & 'C' %in% names(eleFreq)) {
      round(eleFreq['N']/eleFreq['C'],2)
    },
    `O/C` = if ('O' %in% names(eleFreq) & 'C' %in% names(eleFreq)) {
      round(eleFreq['O']/eleFreq['C'],2)
    },
    `P/C` = if ('P' %in% names(eleFreq) & 'C' %in% names(eleFreq)) {
      round(eleFreq['P']/eleFreq['C'],2)
    },
    `S/C` = if ('S' %in% names(eleFreq) & 'C' %in% names(eleFreq)) {
      round(eleFreq['S']/eleFreq['C'],2)
    },
    `N/O` = if ('N' %in% names(eleFreq) & 'O' %in% names(eleFreq)) {
      round(eleFreq['N']/eleFreq['O'],2)
    },
    `P/O` = if ('P' %in% names(eleFreq) & 'O' %in% names(eleFreq)) {
      round(eleFreq['P']/eleFreq['O'],2)
    },
    `S/O` = if ('S' %in% names(eleFreq) & 'O' %in% names(eleFreq)) {
      round(eleFreq['S']/eleFreq['O'],2)
    },
    `O/P` = if ('O' %in% names(eleFreq) & 'P' %in% names(eleFreq)) {
      round(eleFreq['O']/eleFreq['P'],2)
    },
    `S/P` = if ('S' %in% names(eleFreq) & 'P' %in% names(eleFreq)) {
      round(eleFreq['S']/eleFreq['P'],2)
    }
  )
  eleRatios <- as.data.frame(eleRatios[!sapply(eleRatios,is.null)])
  score <- data.frame(MF = mf,eleRatios)
  if ('H.C' %in% colnames(score)) {
    score$H.C <- abs(score$H.C - 1.6)
  }
  if ('O.C' %in% colnames(score)) {
    score$O.C <- abs(score$O.C - 0.3)
  }
  if  (T %in% is.nan(score$N.O)) {
    score$N.O[which(is.nan(score$N.O))] <- 0
  }
  if (T %in% is.infinite(score$N.O)) {
    score$N.O[which(is.infinite(score$N.O))] <- eleFreq['N'][which(is.infinite(score$N.O))]
  }
  if (T %in% is.nan(score$P.O)) {
    score$P.O[which(is.nan(score$P.O))] <- 0
  }
  if (T %in% is.infinite(score$P.O)) {
    score$P.O[which(is.infinite(score$P.O))] <- eleFreq['P'][which(is.infinite(score$P.O))]
  }
  if (T %in% is.nan(score$S.O)) {
    score$S.O[which(is.nan(score$S.O))] <- 0
  }
  if (T %in% is.infinite(score$S.O)) {
    score$S.O[which(is.infinite(score$S.O))] <- eleFreq['S'][which(is.infinite(score$S.O))]
  }
  if (T %in% is.nan(score$O.P)) {
    score$O.P[which(is.nan(score$O.P))] <- 0
  }
  if ('O.P' %in% colnames(score)) {
    score$O.P[which(score$O.P >= 3)] <- 0
  }
  if (T %in% is.nan(score$P.S)) {
    score$P.S[which(is.nan(score$P.S))] <- 0
  }
  if (T %in% is.nan(score$S.P)) {
    score$S.P[which(is.nan(score$S.P))] <- 0
  }
  if (T %in% is.infinite(score$S.P)) {
    score$S.P[which(is.infinite(score$S.P))] <- eleFreq['S'][which(is.infinite(score$S.P))]
  }
  score <- data.frame(score,Score = apply(score[,2:ncol(score)],1,sum))
  return(score$Score[1])
}
