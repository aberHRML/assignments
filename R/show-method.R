#' show-AssignmentParameters
#' @description show method for AssignmentParameters class.
#' @param object S4 object of class AssignmentParameters
#' @importFrom methods show
#' @importFrom crayon yellow
#' @importFrom purrr map
#' @export

setMethod('show',signature = 'AssignmentParameters',
          function(object){
            cat(yellow('\nAssignment Parameters:'),'\n')
            cat('\t','No. of cores:\t\t',object@nCores,'\n')
            cat('\t','Cluster type:\t\t',object@clusterType,'\n')
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

#' show-Assignment
#' @description show mehtod for Assignment class.
#' @param object S4 object of class Assignment
#' @importFrom crayon blue red green
#' @importFrom purrr map_dbl
#' @importFrom utils packageVersion
#' @importFrom igraph E
#' @export

setMethod('show',signature = 'Assignment',
          function(object){
            cat(blue('\nMFassign'),red(str_c('v',packageVersion('MFassign') %>% as.character())),'\n')
            cat(yellow('Assignment:'),'\n')
            cat('\t','Features:\t\t',ncol(object@data),'\n')
            cat('\t','Correlations:\t\t',nrow(object@correlations),'\n')
            cat('\t','Relationships:\t\t',nrow(object@relationships),'\n')
            cat('\n')
            if (length(object@addIsoAssign) > 0) {
              cat('\t',green('Adduct & isotope assignment:'),'\n')
              cat('\t\t','MFs:\t\t',length(unique(object@addIsoAssign$assigned$MF)),'\n')
              cat('\t\t','Relationships:\t',object@addIsoAssign$filteredGraph %>% E() %>% length(),'\n')
              cat('\t\t','Assigned:\t',nrow(object@addIsoAssign$assigned),'\n')
              cat('\n')
            }
            if (length(object@transAssign) > 0) {
              cat('\t',green('Transformation assignment:'),'\n')
              cat('\t\t','Iterations:\t',length(object@transAssign),'\n')
              transAssigned <- object@transAssign %>%
                {.[map_dbl(.,length) > 0]} %>%
                map_dbl(~{
                return(nrow(.$assigned))
              }) %>%
                sum()
              cat('\t\t','Assigned:\t',transAssigned,'\n') 
              cat('\n')
            }
            if (nrow(object@assignments) > 0) {
              cat('\t','Total assignments:\t',blue(nrow(object@assignments)),
                  blue(str_c('(',round(nrow(object@assignments)/ncol(object@data) * 100),'%)')),
                        '\n')
              cat('\t','Unique MFs:\t\t',blue(length(unique(object@assignments$MF))),'\n')
              cat('\n') 
            }
          }
)