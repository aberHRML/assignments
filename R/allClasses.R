#' AssignmentParameters
#' @description An S4 class to store assignment parameters.
#' @slot technique assignment technique to use
#' @slot maxM maximum M for which to assign molecular formulas
#' @slot maxMFscore threshold for molecular formula score
#' @slot ppm ppm threshold
#' #' @slot adducts named list of character vectors containing the adducuts to use for each mode
#' @slot limit amu deviation limit for relationship prediction
#' @slot RTwindow retention time window for chromatographic associations
#' @slot adducts list of character vectors containing the adducts to use. List element names should denote ionisation mode.
#' @slot isotopes character vector of isotopes to use
#' @slot transformations character vector of transformations to use
#' @slot adductRules tibble containing adduct formation rules as returned by mzAnnotation::adducts()
#' @slot isotopeRules tibble containing isotope rules as returned by mzAnnotation::isotopes()
#' @slot transformationRules tibble containing transformation rules as returned by mzAnnotation::transformations()
#' @slot nCores number of cores to use for parallisation
#' @slot clusterType cluster type to use for parallisation
#' @export

setClass('AssignmentParameters',
         slots = list(
           technique = 'character',
           maxM = 'numeric',
           maxMFscore = 'numeric',
           ppm = 'numeric',
           limit = 'numeric',
           RTwindow = 'numeric',
           adducts = 'list',
           isotopes = 'character',
           transformations = 'character',
           adductRules = 'tbl_df',
           isotopeRules = 'tbl_df',
           transformationRules = 'tbl_df',
           nCores = 'numeric',
           clusterType = 'character'
         ))

#' Assignment
#' @description An S4 class to store assignment results
#' @slot log list containing assignment logs
#' @slot flags charactor vector containing completed assignment elements
#' @slot parameters An S4 object of class AssignmentParameters containing the assignment parameters
#' @slot correlations A tibble containing the correlations
#' @slot relationships A tibble containing the predicted relationships
#' @slot addIsoAssign A list containing the results of the adduct and isotope assignment
#' @slot transAssign A list containing the results of the transformation assignment
#' @slot assignments A tibble containing the assigned molecular formulas
#' @export

setClass('Assignment',
         slots = list(
           log = 'list',
           flags = 'character',
           parameters = 'AssignmentParameters',
           correlations = 'tbl_df',
           relationships = 'tbl_df',
           addIsoAssign = 'list',
           transAssign = 'list',
           assignments  = 'tbl_df'
         ))