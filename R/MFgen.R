#' @importFrom mzAnnotation generateMF
#' @importFrom tibble as_tibble

MFgen <- function(M,mz,ppm = 6){
  carb <- round(M/12)
  Hs <- round(carb * 2)
  NO <- round(carb / 2)
  PS <- round(carb / 4)
  
  maxi <- c(C = carb,
            H = Hs,
            N = NO,
            O = NO,
            P = PS,
            S = PS)
  
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
  
  res <- generateMF(M,ppm = ppm,charge = 0,validation = gr,composition = maxi) %>% 
    rename(`Theoretical M` = Mass) %>%
    mutate(`Measured M` = M, `Measured m/z` = mz)
  return(res)
}