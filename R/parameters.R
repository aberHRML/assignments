#' S4 class for assignment parameters
#' @rdname AssignmentParameters-class
#' @description An S4 class to store assignment parameters.
#' @slot technique assignment technique to use
#' @slot correlations_parameters list of correlation parameters to be passed to `metabolyseR::correlation()`
#' @slot max_M maximum M for which to assign molecular formulas
#' @slot MF_rank_threshold rank threshold for molecular formula selection
#' @slot ppm ppm threshold
#' @slot limit amu deviation limit for relationship prediction
#' @slot RT_diff_limit limit for retention time differences for correlated features in adduct and isotopic assignment
#' @slot adducts list of character vectors containing the adducts to use. List element names should denote ionisation mode.
#' @slot isotopes character vector of isotopes to use
#' @slot transformations character vector of transformations to use
#' @slot adduct_rules tibble containing adduct formation rules as returned by `mzAnnotation::adducts()`
#' @slot isotope_rules tibble containing isotope rules as returned by `mzAnnotation::isotopes()`
#' @slot transformation_rules tibble containing transformation rules as returned by `mzAnnotation::transformations()`
#' @importFrom mzAnnotation adduct_rules isotope_rules transformation_rules
#' @export

setClass('AssignmentParameters',
         slots = list(
           technique = 'character',
           correlations_parameters = 'list',
           max_M = 'numeric',
           MF_rank_threshold = 'numeric',
           ppm = 'numeric',
           limit = 'numeric',
           RT_diff_limit = 'numeric',
           adducts = 'list',
           isotopes = 'character',
           transformations = 'character',
           adduct_rules = 'tbl_df',
           isotope_rules = 'tbl_df',
           transformation_rules = 'tbl_df'
         ),
         prototype = list(
           technique = 'FIE-HRMS',
           correlations_parameters = list(method = 'spearman',
                                          pAdjustMethod = 'bonferroni',
                                          corPvalue = 0.05,
                                          minCoef = 0.7,
                                          maxCor = 500000),
           max_M = 800,
           MF_rank_threshold = 3,
           ppm = 6,
           limit = 0.001,
           RT_diff_limit = numeric(),
           isotopes = c('13C','18O','13C2'),
           adducts = list(n = c("[M-H]1-", "[M+Cl]1-", "[M+K-2H]1-", 
                                "[M-2H]2-", "[M+Cl37]1-","[2M-H]1-"),
                          p = c('[M+H]1+','[M+K]1+','[M+Na]1+','[M+K41]1+',
                                '[M+2H]2+','[2M+H]1+')),
           transformations = transformation_rules()$`MF Change`,
           adduct_rules = adduct_rules(),
           isotope_rules = isotope_rules(),
           transformation_rules = transformation_rules()
         ))

setValidity('AssignmentParameters',function(object){
  technique_correct <- technique(object) %in% availableTechniques()
  
  if (isFALSE(technique_correct)) {
    availableTechniques() %>% 
      paste(collapse = ', ') %>% 
      paste0('Technique should be one of ',.)
  }
  else TRUE
})

#' @importFrom metabolyseR correlationsParameters

setValidity('AssignmentParameters',function(object){
  
  correlations_parameters <- correlationsParameters() %>% 
    names()
  
  if (!any(names(object@correlations_parameters) %in% 
           correlations_parameters)) {
    correlations_parameters %>% 
      paste0('`',.,'`') %>% 
      paste(collapse = ', ') %>% 
      paste0('Correlations parameters should only include ',.)
  }
  else TRUE
})

#' @importFrom methods show
#' @importFrom crayon yellow
#' @importFrom purrr map

setMethod('show',signature = 'AssignmentParameters',
          function(object){
            cat(yellow('\nAssignment Parameters:'),'\n')
            cat('\n')
            cat('\t','Technique:\t\t',object@technique,'\n')
            cat('\t','Max M:\t\t\t',object@max_M,'\n')
            cat('\t','MF rank threshold:\t',object@MF_rank_threshold,'\n')
            cat('\t','PPM threshold:\t\t',object@ppm,'\n')
            cat('\t','Relationship limit:\t',object@limit,'\n')
            
            if (object@technique != 'FIE') {
              cat('\t','RT limit:\t\t',object@RT_diff_limit,'\n')
            }
            
            
            cat('\t','Correlations:\n')
            correlations(object) %>% 
              paste0('\t\t',names(.),': ',.,'\n') %>% 
              cat()
            
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
#' ## Return max M
#' maxM(assignment_parameters)
#' 
#' ## Set max M
#' maxM(assignment_parameters) <- 500
#' 
#' ## Return MF rank threshold
#' MFrankThreshold(assignment_parameters)
#' 
#' ## Set MF rank threshold
#' MFrankThreshold(assignment_parameters) <- 3
#' 
#' ## Return ppm
#' ppm(assignment_parameters)
#' 
#' ## Set ppm
#' ppm(assignment_parameters) <- 3
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
#' adductRules(assignment_parameters) <- mzAnnotation::adduct_rules()
#' 
#' ## Return isotope rules
#' isotopeRules(assignment_parameters)
#' 
#' ## Set isotope rules
#' isotopeRules(assignment_parameters) <- mzAnnotation::isotope_rules()
#' 
#' ## Return transformation rules
#' transformationRules(assignment_parameters)
#' 
#' ## Set transformation rules
#' transformationRules(assignment_parameters) <- mzAnnotation::transformation_rules()
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

setGeneric('correlations',
           function(x) standardGeneric('correlations'))

#' @rdname parameters

setMethod('correlations',signature = 'AssignmentParameters',
          function(x) x@correlations_parameters)

#' @rdname parameters
#' @export

setGeneric('correlations<-',
           function(x,value) standardGeneric('correlations<-'))

#' @rdname parameters
#' @importFrom methods validObject

setMethod('correlations<-',signature = c('AssignmentParameters','list'),
          function(x,value){
            x@correlations_parameters <- value
            validObject(x)
            return(x)
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

setGeneric('maxM',function(x)
  standardGeneric('maxM'))

#' @rdname parameters

setMethod('maxM',signature = 'AssignmentParameters',
          function(x){
            x@max_M
          })

#' @rdname parameters
#' @export

setGeneric('maxM<-',function(x,value)
  standardGeneric('maxM<-'))

#' @rdname parameters

setMethod('maxM<-',signature = 'AssignmentParameters',
          function(x,value){
            x@max_M <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('MFrankThreshold',function(x)
  standardGeneric('MFrankThreshold'))

#' @rdname parameters

setMethod('MFrankThreshold',signature = 'AssignmentParameters',
          function(x){
            x@MF_rank_threshold
          })

#' @rdname parameters
#' @export

setGeneric('MFrankThreshold<-',function(x,value)
  standardGeneric('MFrankThreshold<-'))

#' @rdname parameters

setMethod('MFrankThreshold<-',signature = 'AssignmentParameters',
          function(x,value){
            x@MF_rank_threshold <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('ppm',function(x)
  standardGeneric('ppm'))

#' @rdname parameters

setMethod('ppm',signature = 'AssignmentParameters',
          function(x){
            x@ppm
          })

#' @rdname parameters
#' @export

setGeneric('ppm<-',function(x,value)
  standardGeneric('ppm<-'))

#' @rdname parameters

setMethod('ppm<-',signature = 'AssignmentParameters',
          function(x,value){
            x@ppm <- value
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
            x@adduct_rules
          })

#' @rdname parameters
#' @export

setGeneric('adductRules<-',function(x,value)
  standardGeneric('adductRules<-'))

#' @rdname parameters

setMethod('adductRules<-',signature = 'AssignmentParameters',
          function(x,value){
            x@adduct_rules <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('isotopeRules',function(x)
  standardGeneric('isotopeRules'))

#' @rdname parameters

setMethod('isotopeRules',signature = 'AssignmentParameters',
          function(x){
            x@isotope_rules
          })

#' @rdname parameters
#' @export

setGeneric('isotopeRules<-',function(x,value)
  standardGeneric('isotopeRules<-'))

#' @rdname parameters

setMethod('isotopeRules<-',signature = 'AssignmentParameters',
          function(x,value){
            x@isotope_rules <- value
            return(x)
          })

#' @rdname parameters
#' @export

setGeneric('transformationRules',function(x)
  standardGeneric('transformationRules'))

#' @rdname parameters

setMethod('transformationRules',signature = 'AssignmentParameters',
          function(x){
            x@transformation_rules
          })

#' @rdname parameters
#' @export

setGeneric('transformationRules<-',function(x,value)
  standardGeneric('transformationRules<-'))

#' @rdname parameters

setMethod('transformationRules<-',signature = 'AssignmentParameters',
          function(x,value){
            x@transformation_rules <- value
            return(x)
          })

#' Available techniques
#' @description Available techniques for molecular formula assignment.
#' @return A `character` vector of technique names.
#' @examples 
#' availableTechniques()
#' @export

availableTechniques <- function(){
  c('FIE-HRMS','RP-LC-HRMS','NP-LC-HRMS')
}

#' Assignment parameters
#' @description Return default assignment parameters for a given technique.
#' @param technique technique to use for assignment
#' @importFrom methods new
#' @export

assignmentParameters <- function(technique = availableTechniques()){
  
  technique <- match.arg(technique,
                         choices = availableTechniques())
  
  parameters <- switch(technique,
                       `FIE-HRMS` = new('AssignmentParameters'),
                       `RP-LC-HRMS` = new('AssignmentParameters',
                                          technique = 'RP-LC-HRMS',
                                          max_M = 700,
                                          RT_diff_limit = 1/60),
                       `NP-LC-HRMS` = new('AssignmentParameters',
                                          technique = 'NP-LC-HRMS',
                                          max_M = 700,
                                          RT_diff_limit = 1/60,
                                          adducts = list(n = c("[M-H]1-", "[M+Cl]1-", "[M+K-2H]1-", 
                                                               "[M-2H]2-", "[M+Cl37]1-","[2M-H]1-"),
                                                         p = c('[M+H]1+','[M+K]1+','[M+Na]1+','[M+K41]1+',
                                                               '[M+NH4]1+','[M+2H]2+','[2M+H]1+'))))
  
  return(parameters)
} 