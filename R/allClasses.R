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
           log = list(date = date(),verbose = TRUE),
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