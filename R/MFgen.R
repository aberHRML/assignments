#' @importFrom mzAnnotation generateMF

MFgen <- function(M,mz,ppm = 6){
  carb <- round(M/12)
  Hs <- round(carb * 2)
  NO <- round(carb / 2)
  PS <- round(carb / 4)
  
  maxi <- c(C = carb,
            iC = 0,
            H = Hs,
            iH = 0,
            N = NO,
            iN = 0,
            O = NO,
            iO = 0,
            F = 0 ,
            Na = 0,
            Si = 0,
            P = PS,
            S = PS,
            Cl = 0,
            iCl = 0,
            Br = 0,
            iBr = 0,
            K = 0,
            iK = 0)
  
  if (M < 100) {
    ppm <- 10
  } else {
    ppm <- ppm
  }
  
  ppm <- (ppm/10^6 * mz)/M * 10^6
  
  if (M < 200) {
    gr <- F
  } else {
    gr <- T
  }
  
  res <- generateMF(M,ppm = ppm,charge = 0,applygr = gr,composition = maxi)
  res$`m/z` <- round(as.numeric(res$`m/z`),5) 
  colnames(res)[2] <- 'Theoretical M'
  res <- mutate(res, 'Measured M' = M, `Measured m/z` = mz)
  return(res)
}