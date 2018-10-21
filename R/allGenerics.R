setGeneric('prepCorrelations', function(assignment){
  standardGeneric('prepCorrelations')
})

setGeneric("relationships", function(assignment) {
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
setGeneric('plotNetwork',function(assignment, layout = 'kk', rThreshold = 0.7, labels = F){
standardGeneric('plotNetwork')
})