#' @export
setClass('AnnotationParameters',
         slots = list(
           technique = 'character',
           ppm = 'numeric',
           limit = 'numeric',
           adducts = 'list',
           nCores = 'numeric'
         ))

#' @export
setClass('Annotation',
         slots = list(
           parameters = 'AnnotationParameters',
           correlations = 'tbl_df',
           relationships = 'tbl_df',
           annotations = 'tbl_df',
           results = 'tbl_df'
         ))