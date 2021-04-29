#' AssignmentParameters
#' @description An S4 class to store assignment parameters.
#' @slot technique assignment technique to use
#' @slot correlations list of correlation parameters to be passed to metabolyseR correlation analysis
#' @slot filter list of r and n thresholds for filtering correlations
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
#' @importFrom mzAnnotation transformations
#' @export

setClass('AssignmentParameters',
         slots = list(
           technique = 'character',
           correlations = 'list',
           filter = 'list',
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
           transformationRules = 'tbl_df'
         ),
         prototype = list(
           technique = 'FIE',
           correlations = list(method = 'pearson',pAdjustMethod = 'bonferroni',corPvalue = 0.05),
           filter = list(rthresh = 0.7,n = 100000,rIncrement = 0.01,nIncrement = 20000),
           maxM = 1000,
           maxMFscore = 5,
           ppm = 5,
           limit = 0.001,
           RTwindow = numeric(),
           isotopes = c('13C','18O','13C2'),
           adducts = list(n = c("[M-H]1-", "[M+Cl]1-", "[M+K-2H]1-", 
                                "[M-2H]2-", "[M+Cl37]1-","[2M-H]1-"),
                          p = c('[M+H]1+','[M+K]1+','[M+Na]1+','[M+K41]1+',
                                '[M+NH4]1+','[M+2H]2+','[2M+H]1+')),
           transformations = transformations()$`MF Change`,
           adductRules = adducts(),
           isotopeRules = isotopes(),
           transformationRules = transformations()
         ))

#' Assignment
#' @description An S4 class to store assignment results
#' @slot log list containing assignment logs
#' @slot flags charactor vector containing completed assignment elements
#' @slot parameters An S4 object of class AssignmentParameters containing the assignment parameters
#' @slot data A tibble containing the peak intensity matrix
#' @slot correlations A tibble containing the correlations
#' @slot preparedCorrelations A tibble containing the prepared correlations ready for analysis
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
           data = 'tbl_df',
           correlations = 'tbl_df',
           preparedCorrelations = 'tbl_df',
           relationships = 'tbl_df',
           addIsoAssign = 'list',
           transAssign = 'list',
           assignments  = 'tbl_df'
         ),
         prototype = list(
           log = list(date = date(),verbose = verbose),
           flags = character(),
           parameters = new('AssignmentParameters'),
           data = tibble(),
           correlations = tibble(),
           preparedCorrelations = tibble(),
           relationships = tibble(),
           addIsoAssign = list(),
           transAssign = list(),
           assignments  = tibble()
         ))