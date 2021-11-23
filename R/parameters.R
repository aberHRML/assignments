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
#' @importFrom mzAnnotation adduct_rules isotope_rules transformation_rules
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
           correlations_parameters = list(method = 'pearson',
                                          pAdjustMethod = 'bonferroni',
                                          corPvalue = 0.05),
           filter = list(rthresh = 0.7,
                         n = 200000,
                         rIncrement = 0.01,
                         nIncrement = 20000),
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
           transformations = transformation_rules()$`MF Change`,
           adductRules = adduct_rules(),
           isotopeRules = isotope_rules(),
           transformationRules = transformation_rules()
         ))

# setValidity('AssignmentParameters',function(object){
#   adducts_present <- adducts(object) %in% adductRules(object)$Name
#   
#   if (FALSE %in% adducts_present){
#     missing_adducts <- adducts(object)[adducts_present] %>% 
#       str_c(collapse = ', ')
#     
#     str_c('Specified adducts ',missing_adducts,' not present in adduct rules.')
#   } else {
#     TRUE
#   }
# })

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
#' ## Return technique
#' technique(assignment_parameters)
#' 
#' ## Return limit
#' limit(assignment_parameters)
#' 
#' ## Set limit
#' limit(assignment_parameters) <- 0.002
#' 
#' ## Return isotopes
#' isotopes(assignment_parameters)
#' 
#' ## Set isotopes
#' isotopes(assignment_parameters) <- '13C'
#' 
#' ## Return adducts
#' adducts(assignment_parameters)
#' 
#' ## Set adducts
#' adducts(assignment_parameters) <- list(n = c('[M-H]1-','[M+Cl]1-'),
#'                                    p = c('[M+H]1+','[M+K]1+'))
#'                                    
#' ## Return transformations
#' transformations(assignment_parameters)
#' 
#' ## Set transformations
#' transformations(assignment_parameters) <- "M - [O] + [NH2]"
#' 
#' ## Return adduct rules
#' adductRules(assignment_parameters)
#' 
#' ## Set adduct rules
#' adductRules(assignment_parameters) <- mzAnnotation::adduct_rules())
#' 
#' ## Return isotope rules
#' isotopeRules(assignment_parameters)
#' 
#' ## Set isotope rules
#' isotopeRules(assignment_parameters) <- mzAnnotation::isotope_rules())
#' 
#' ## Return transformation rules
#' transformationRules(assignment_parameters)
#' 
#' ## Set transformation rules
#' transformationRules(assignment_parameters) <- mzAnnotation::transformation_rules())
#' @export

setGeneric('technique',function(x)
  standardGeneric('technique'))

#' @rdname parameters

setMethod('technique',signature = 'AssignmentParameters',
          function(x){
            x@technique
          })

#' @rdname parameters
#' @export

setGeneric('limit',function(x)
  standardGeneric('limit'))

#' @rdname parameters

setMethod('limit',signature = 'AssignmentParameters',
          function(x){
            x@limit
          })

#' @rdname parameters
#' @export

setGeneric('limit<-',function(x,value)
  standardGeneric('limit<-'))

#' @rdname parameters

setMethod('limit<-',signature = 'AssignmentParameters',
          function(x,value){
            x@limit <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('isotopes',function(x)
  standardGeneric('isotopes'))

#' @rdname parameters

setMethod('isotopes',signature = 'AssignmentParameters',
          function(x){
            x@isotopes
          })

#' @rdname parameters
#' @export

setGeneric('isotopes<-',function(x,value)
  standardGeneric('isotopes<-'))

#' @rdname parameters

setMethod('isotopes<-',signature = 'AssignmentParameters',
          function(x,value){
            x@isotopes <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('adducts',function(x)
  standardGeneric('adducts'))

#' @rdname parameters

setMethod('adducts',signature = 'AssignmentParameters',
          function(x){
            x@adducts
          })

#' @rdname parameters
#' @export

setGeneric('adducts<-',function(x,value)
  standardGeneric('adducts<-'))

#' @rdname parameters

setMethod('adducts<-',signature = 'AssignmentParameters',
          function(x,value){
            x@adducts <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('transformations',function(x)
  standardGeneric('transformations'))

#' @rdname parameters

setMethod('transformations',signature = 'AssignmentParameters',
          function(x){
            x@transformations
          })

#' @rdname parameters
#' @export

setGeneric('transformations<-',function(x,value)
  standardGeneric('transformations<-'))

#' @rdname parameters

setMethod('transformations<-',signature = 'AssignmentParameters',
          function(x,value){
            x@transformations <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('adductRules',function(x)
  standardGeneric('adductRules'))

#' @rdname parameters

setMethod('adductRules',signature = 'AssignmentParameters',
          function(x){
            x@adductRules
          })

#' @rdname parameters
#' @export

setGeneric('adductRules<-',function(x,value)
  standardGeneric('adductRules<-'))

#' @rdname parameters

setMethod('adductRules<-',signature = 'AssignmentParameters',
          function(x,value){
            x@adductRules <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('isotopeRules',function(x)
  standardGeneric('isotopeRules'))

#' @rdname parameters

setMethod('isotopeRules',signature = 'AssignmentParameters',
          function(x){
            x@isotopeRules
          })

#' @rdname parameters
#' @export

setGeneric('isotopeRules<-',function(x,value)
  standardGeneric('isotopeRules<-'))

#' @rdname parameters

setMethod('isotopeRules<-',signature = 'AssignmentParameters',
          function(x,value){
            x@isotopeRules <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('transformationRules',function(x)
  standardGeneric('transformationRules'))

#' @rdname parameters

setMethod('transformationRules',signature = 'AssignmentParameters',
          function(x){
            x@isotopeRules
          })

#' @rdname parameters
#' @export

setGeneric('transformationRules<-',function(x,value)
  standardGeneric('transformationRules<-'))

#' @rdname parameters

setMethod('transformationRules<-',signature = 'AssignmentParameters',
          function(x,value){
            x@transformationRules <- value
            return(x)
          })