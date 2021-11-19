setGeneric('calcCorrelations', function(assignment){
  standardGeneric('calcCorrelations')
})

setGeneric('filterCorrelations', function(assignment){
  standardGeneric('filterCorrelations')
})

setGeneric('prepCorrelations', function(assignment){
  standardGeneric('prepCorrelations')
})

setGeneric("relationships", function(assignment,transformations = TRUE) {
  standardGeneric("relationships")
})

setGeneric("addIsoAssign", function(assignment) {
  standardGeneric("addIsoAssign")
})

setGeneric("transformationAssign", function(assignment) {
  standardGeneric("transformationAssign")
})

#' @rdname summariseAssignment
setGeneric('summariseAssignment',function(assignment){
  standardGeneric('summariseAssignment')
})

#' @rdname plotSpectrum
setGeneric('plotSpectrum',function(assignment,MF){
  standardGeneric('plotSpectrum')
})