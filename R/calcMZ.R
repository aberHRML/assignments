
calcMZ <- function(M,isotope,adduct){
  
  if (!is.na(isotope)) {
    isorule <- rules$isotopes[which(rules$isotopes == isotope),]
    M <- M + isorule$Mass.Difference
  }
  
  addrule <- rules$adducts[which(rules$adducts$Name == adduct),]

  mz <- (M * addrule$xM) / addrule$Charge + addrule$Add
  
  return(round(mz,5))
}