
calcM <- function(mz,isotope,adduct){

    addrule <- rules$adducts[which(rules$adducts$Name == adduct),]
    
    M <- ((mz - addrule$Add) * addrule$Charge) / addrule$xM 
    
    if (!is.na(isotope)) {
      isorule <- rules$isotopes[which(rules$isotopes == isotope),]
      M <- M - isorule$Mass.Difference
    }
    return(round(M,5))
}