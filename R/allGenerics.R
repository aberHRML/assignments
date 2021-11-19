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

#' @rdname assignments
setGeneric('assignments',function(assignment){
  standardGeneric('assignments')
})

#' @rdname summariseAssignment
setGeneric('summariseAssignment',function(assignment){
  standardGeneric('summariseAssignment')
})

#' @rdname assignmentData
setGeneric('assignmentData',function(assignment){
  standardGeneric('assignmentData')
})

#' @rdname assignedData
setGeneric('assignedData',function(assignment){
  standardGeneric('assignedData')
})

#' @rdname plotAdductDist
setGeneric('plotAdductDist',function(assignment){
  standardGeneric('plotAdductDist')
})

#' @rdname plotSpectrum
setGeneric('plotSpectrum',function(assignment,MF){
  standardGeneric('plotSpectrum')
})