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

setGeneric('doAssignment',function(assignment){
  standardGeneric('doAssignment')
})

#' @rdname assignments
setGeneric('assignments',function(assignment){
  standardGeneric('assignments')
})

#' @rdname summariseAssignment
setGeneric('summariseAssignment',function(assignment){
  standardGeneric('summariseAssignment')
})

#' @rdname plotNetwork
setGeneric('plotNetwork',function(assignment, layout = 'stress', rThreshold = 0.7){
standardGeneric('plotNetwork')
})

#' @rdname plotAdductDist
setGeneric('plotAdductDist',function(assignment){
  standardGeneric('plotAdductDist')
})

#' @rdname plotFeatureSolutions
setGeneric('plotFeatureSolutions',function(assignment,feature,maxComponents = 10){
  standardGeneric('plotFeatureSolutions')
})

#' @rdname plotSpectrum
setGeneric('plotSpectrum',function(assignment,MF){
  standardGeneric('plotSpectrum')
})