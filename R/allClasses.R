#' @export
setClass('AssignmentParameters',
         slots = list(
           technique = 'character',
           maxM = 'numeric',
           maxMFscore = 'numeric',
           ppm = 'numeric',
           limit = 'numeric',
           isotopes = 'character',
           adducts = 'list',
           nCores = 'numeric',
           clusterType = 'character'
         ))

#' @export
setClass('Assignment',
         slots = list(
           parameters = 'AssignmentParameters',
           correlations = 'tbl_df',
           relationships = 'tbl_df',
           addIsoAssign = 'list',
           transAssign = 'list',
           assignments  = 'tbl_df'
         ))