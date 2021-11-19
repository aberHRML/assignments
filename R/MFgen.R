maxElementFrequencies <- function(M){
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
  
  return(maxi)
}


#' @importFrom mzAnnotation generateMF
#' @importFrom tibble as_tibble

MFgen <- function(M,mz,ppm = 6){
  
  max_element_frequencies <- maxElementFrequencies(M)
  
  if (M < 100) {
    ppm <- 10
  } else {
    ppm <- ppm
  }
  
  ppm <- (ppm/10^6 * mz)/M * 10^6
  
  if (M < 200) {
    gr <- FALSE
  } else {
    gr <- TRUE
  }
  
  res <- generateMF(M,ppm = ppm,
                    charge = 0,
                    validation = gr,
                    element_max = max_element_frequencies) %>% 
    rename(`Theoretical M` = Mass) %>%
    mutate(`Measured M` = M, `Measured m/z` = mz)
  return(res)
}