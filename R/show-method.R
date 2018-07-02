#' show-AssignmentParameters
#' @description show method for AssignmentParameters class.
#' @param object S4 object of class AssignmentParameters
#' @importFrom methods show
#' @importFrom crayon yellow
#' @export

setMethod('show',signature = 'AssignmentParameters',
          function(object){
            cat(yellow('\nAssignment Parameters:'),'\n')
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