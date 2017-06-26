#' @export
setClass('AnnotationParameters',
         slots = list(
           technique = 'character',
           ppm = 'numeric',
           nCores = 'numeric'
         ))

#' @export
setClass('Annotation',
         slots = list(
           correlations = 'tbl_df',
           relationships = 'tbl_df',
           annotations = 'tbl_df',
           results = 'list'
         ))