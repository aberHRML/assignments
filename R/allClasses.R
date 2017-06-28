#' @export
setClass('AnnotationParameters',
         slots = list(
           technique = 'character',
           maxM = 'numeric',
           maxMFscore = 'numeric',
           ppm = 'numeric',
           limit = 'numeric',
           isotopes = 'character',
           adducts = 'list',
           nCores = 'numeric'
         ))

#' @export
setClass('Annotation',
         slots = list(
           parameters = 'AnnotationParameters',
           correlations = 'tbl_df',
           relationships = 'tbl_df',
           addIsoAssign = 'list',
           transAssign = 'list',
           assignments  = 'tbl_df'
         ))