#' show-Assignment
#' @description show mehtod for Assignment class.
#' @param object S4 object of class Assignment
#' @importFrom crayon blue red green
#' @importFrom purrr map_dbl
#' @importFrom utils packageVersion
#' @importFrom igraph E
#' @exporti

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