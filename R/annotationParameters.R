#' annotationParameters
#' @importFrom parallel detectCores
#' @export

annotationParameters <- function(technique = NULL){
  if (is.null(technique)) {
    cat('\nAvailable Techniques:','\t FIE')
  }
  
  if (technique == 'FIE') {
    new('AnnotationParameters',
        technique = 'FIE',
        ppm = 5,
        nCores = detectCores()
        )
  }
} 