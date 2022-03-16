#' Assignment
#' @rdname Assignment
#' @description An S4 class to store assignment results
#' @slot log list containing assignment logs
#' @slot flags charactor vector containing completed assignment elements
#' @slot data A tibble containing the peak intensity matrix
#' @slot correlations A tibble containing the correlations
#' @slot preparedCorrelations A tibble containing the prepared correlations ready for analysis
#' @slot relationships A tibble containing the predicted relationships
#' @slot addIsoAssign A list containing the results of the adduct and isotope assignment
#' @slot transAssign A list containing the results of the transformation assignment
#' @slot assignments A tibble containing the assigned molecular formulas
#' @export

setClass('Assignment',
         contains = 'AssignmentParameters',
         slots = list(
           log = 'list',
           flags = 'character',
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
           data = tibble(),
           correlations = tibble(),
           preparedCorrelations = tibble(),
           relationships = tibble(),
           addIsoAssign = list(),
           transAssign = list(),
           assignments  = tibble()
         ))

#' @importFrom crayon blue red green
#' @importFrom purrr map_dbl
#' @importFrom utils packageVersion
#' @importFrom igraph E

setMethod('show',signature = 'Assignment',
          function(object){
            cat(blue('\nMFassign'),red(str_c('v',packageVersion('MFassign') %>% as.character())),'\n')
            cat(yellow('Assignment:'),'\n')
            cat('\t','Features:\t\t',ncol(object@data),'\n')
            cat('\t','Correlations:\t\t',nrow(object@preparedCorrelations),'\n')
            cat('\t','Relationships:\t\t',nrow(relationships(object)),'\n')
            cat('\n')
            if (length(object@addIsoAssign) > 0) {
              cat('\t',green('Adduct & isotope assignment:'),'\n')
              cat('\t\t','MFs:\t\t',length(unique(object@addIsoAssign$assigned$MF)),'\n')
              cat('\t\t','Relationships:\t',object@addIsoAssign$filtered_graph %>% E() %>% length(),'\n')
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

#' Assignment accessors
#' @rdname accessors
#' @description Access methods for Assignment S4 class
#' @param assignment S4 object of class Assignment
#' @export

setGeneric('relationships',function(assignment)
  standardGeneric('relationships'))

#' @rdname accessors

setMethod('relationships',signature = 'Assignment',
          function(assignment){
            assignment@relationships
          })

setGeneric('relationships<-',function(assignment,value)
  standardGeneric('relationships<-'))

setMethod('relationships<-',signature = 'Assignment',
          function(assignment,value){
            assignment@relationships <- value
            return(assignment)
          })

#' @rdname accessors
#' @export

setGeneric('assignments',function(assignment)
  standardGeneric('assignments'))

#' @rdname accessors

setMethod('assignments',signature = 'Assignment',
          function(assignment){
            assignment@assignments
          })

#' assignmentData
#' @rdname assignmentData
#' @description Return data table used for assignments.
#' @param assignment S4 object of class Assignment
#' @export

setGeneric('assignmentData',function(assignment)
  standardGeneric('assignmentData'))

#' @rdname assignmentData

setMethod('assignmentData', signature = 'Assignment',
          function(assignment){
            assignment@data
          })

#' assignedData
#' @rdname assignedData
#' @description Return data table used for assignments with feature assignments added to column names.
#' @param assignment S4 object of class Assignment
#' @export

setGeneric('assignedData',function(assignment)
  standardGeneric('assignedData'))

#' @rdname assignedData

setMethod('assignedData', signature = 'Assignment',
          function(assignment){
            
            d <- assignment %>%
              assignmentData()
            
            assignedFeats <- assignment %>%
              assignments() %>%
              select(Feature,Name)
            
            assignedFeats <- left_join(
              tibble(Feature = d %>% colnames()),
              assignedFeats, 
              by = "Feature")
            
            assignedFeats$Name[is.na(assignedFeats$Name)] <- assignedFeats$Feature[is.na(assignedFeats$Name)] 
            
            assignedFeats <- assignedFeats %>%
              filter(!duplicated(Feature))
            
            colnames(d) <- assignedFeats$Name 
            
            return(d)
          })

#' Summarise assignments
#' @rdname summariseAssignment
#' @description Summarise features assigned to molecular formulas.
#' @param assignment S4 object of class Assignment
#' @importFrom dplyr desc
#' @export

setGeneric('summariseAssignment',function(assignment)
  standardGeneric('summariseAssignment'))

#' @rdname summariseAssignment

setMethod('summariseAssignment',signature = 'Assignment',
          function(assignment){
            assigned <- assignment %>%
              assignments() %>%
              split(.$MF) %>%
              map(~{
                d <- .
                d$Isotope[is.na(d$Isotope)] <- ''
                d <- d %>%
                  mutate(IIP = str_c(Isotope,Adduct,sep = ' ')) %>%
                  arrange(`Measured m/z`)
                tibble(MF = d$MF[1],Features = str_c(d$Feature,collapse = '; '),`Isotopes & Ionisation Products` = str_c(d$IIP,collapse = '; '),Count = nrow(d))
              }) %>%
              bind_rows() %>%
              arrange(desc(Count))
            return(assigned)
          })