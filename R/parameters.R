#' S4 class for assignment parameters
#' @rdname AssignmentParameters-class
#' @description An S4 class to store assignment parameters.
#' @slot technique assignment technique to use
#' @slot correlations_parameters list of correlation parameters to be passed to metabolyseR correlation analysis
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
           correlations_parameters = 'list',
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
           correlations_parameters = list(method = 'pearson',pAdjustMethod = 'bonferroni',corPvalue = 0.05),
           filter = list(rthresh = 0.7,n = 200000,rIncrement = 0.01,nIncrement = 20000),
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

#' @importFrom methods show
#' @importFrom crayon yellow
#' @importFrom purrr map

setMethod('show',signature = 'AssignmentParameters',
          function(object){
            cat(yellow('\nAssignment Parameters:'),'\n')
            cat('\n')
            cat('\t','Technique:\t\t',object@technique,'\n')
            cat('\t','Max M:\t\t\t',object@maxM,'\n')
            cat('\t','Max MF score:\t\t',object@maxMFscore,'\n')
            cat('\t','PPM threshold:\t\t',object@ppm,'\n')
            cat('\t','Relationship limit:\t',object@limit,'\n')
            
            if (object@technique != 'FIE') {
              cat('\t','RT window:\t\t',object@RTwindow,'\n')
            }
            
            cat('\n\t','Adducts:','\n')
            adducts <- map(names(object@adducts),~{
              a <- str_c(object@adducts[[.]],collapse = ', ')
              str_c(.,': ',a)
            }) %>%
              str_c(collapse = '\n\t ')
            cat('\t',adducts,'\n')
            
            cat('\t','Isotopes:',str_c(object@isotopes,collapse = ', '),'\n')
            
            cat('\t','Transformations:',str_c(object@transformations,collapse = ', '))
            
            cat('\n')
          }
)

#' Parameter get and set methods
#' @rdname parameters
#' @description Get and set methods for the `AssignmentParameters` S4 class.
#' @param x S4 object of class `AssignmentParameters`
#' @param value value to set
#' @examples 
#' assignment_parameters <- assignmentParameters('FIE')
#' 
#' ## Return isotopes
#' iso(assignment_parameters)
#' 
#' ## Set isotopes
#' iso(assignment_parameters) <- '13C'
#' 
#' ## Return adducts
#' add(assignment_parameters)
#' 
#' ## Set adducts
#' add(assignment_parameters) <- list(n = c('[M-H]1-','[M+Cl]1-'),
#'                                    p = c('[M+H]1+','[M+K]1+'))
#' 
#' @export

setGeneric('iso',function(x)
  standardGeneric('iso'))

#' @rdname parameters

setMethod('iso',signature = 'AssignmentParameters',
          function(x){
            x@isotopes
          })

#' @rdname parameters
#' @export

setGeneric('iso<-',function(x,value)
  standardGeneric('iso<-'))

#' @rdname parameters

setMethod('iso<-',signature = 'AssignmentParameters',
          function(x,value){
            x@isotopes <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('add',function(x)
  standardGeneric('add'))

#' @rdname parameters

setMethod('add',signature = 'AssignmentParameters',
          function(x){
            x@adducts
          })

#' @rdname parameters
#' @export

setGeneric('add<-',function(x,value)
  standardGeneric('add<-'))

#' @rdname parameters

setMethod('add<-',signature = 'AssignmentParameters',
          function(x,value){
            x@adducts <- value
            return(x)
          })
