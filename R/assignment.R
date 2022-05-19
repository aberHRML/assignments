#' Assignment
#' @rdname Assignment-class
#' @description An S4 class to store assignment results
#' @slot log list containing assignment logs
#' @slot flags character vector containing completed assignment elements
#' @slot data A tibble containing the peak intensity matrix
#' @slot correlations A tibble containing the correlations
#' @slot relationships A tibble containing the predicted relationships
#' @slot addIsoAssign A list containing the results of the adduct and isotope assignment
#' @slot transAssign A list containing the results of the transformation assignment
#' @slot assignments A tibble containing the assigned molecular formulas

setClass('Assignment',
         contains = 'AssignmentParameters',
         slots = list(
           log = 'list',
           data = 'tbl_df',
           correlations = 'tbl_df',
           relationships = 'tbl_df',
           addIsoAssign = 'list',
           transAssign = 'list',
           assignments  = 'tbl_df'
         ),
         prototype = list(
           log = list(date = date(),verbose = TRUE),
           data = tibble(),
           correlations = tibble(),
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
            cat(blue('\nassignments'),red(str_c('v',packageVersion('assignments') %>% as.character())),'\n')
            cat(yellow('Assignment:'),'\n')
            cat('\t','Features:\t\t',ncol(object@data),'\n')
            cat('\t','Correlations:\t\t',nrow(object@correlations),'\n')
            cat('\t','Relationships:\t\t',nrow(relationships(object)),'\n')
            cat('\n')
            if (length(object@addIsoAssign) > 0) {
              addIsoAssigned <- object %>% 
                assignments() %>% 
                filter(str_detect(Iteration,'A&I')) 
              cat('\t',green('Adduct & isotope assignment:'),'\n')
              cat('\t\t','Iterations:\t',length(object@addIsoAssign),'\n')
              cat('\t\t',
                  'MFs:\t\t',
                  addIsoAssigned %>% 
                    select(MF) %>% 
                    distinct() %>% 
                    nrow(),
                  '\n')
              cat('\t\t','Assigned:\t',nrow(addIsoAssigned),'\n') 
              cat('\n')
            }
            if (length(object@transAssign) > 0) {
              transAssign <- object %>% 
                assignments() %>% 
                filter(str_detect(Iteration,'T')) 
              cat('\t',green('Transformation assignment:'),'\n')
              cat('\t\t','Iterations:\t',length(object@transAssign),'\n')
              cat('\t\t',
                  'MFs:\t\t',
                  transAssign %>% 
                    select(MF) %>% 
                    distinct() %>% 
                    nrow(),
                  '\n')
              cat('\t\t','Assigned:\t',nrow(transAssign),'\n') 
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
#' @param iteration the assignment iteration
#' @param type the graph type to return. `filtered` returns the assignment graph after component selection. `all` returns all assignment components.
#' @examples 
#' assignment <- new('Assignment',
#'                   data = feature_data)
#' 
#' ## Return feature data
#' featureData(assignment)
#' 
#' ## Return correlations
#' correlations(assignment)
#'
#' ## Return relationships
#' relationships(assignment)
#' 
#' ## Return the available iterations
#' iterations(assignment)
#'
#' ## Return a selected graph
#' \dontrun{
#' graph(assignment,'A&I1')
#' }
#' 
#' ## Return assignments
#' assignments(assignment)
#' @export

setGeneric('featureData',function(assignment)
  standardGeneric('featureData'))

#' @rdname accessors

setMethod('featureData', signature = 'Assignment',
          function(assignment){
            assignment@data
          })

#' @rdname accessors
#' @export

setGeneric('correlations',function(assignment)
  standardGeneric('correlations'))

#' @rdname accessors

setMethod('correlations',signature = 'Assignment',
          function(assignment){
            assignment@relationships
          })

#' @rdname accessors
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

setGeneric('iterations',function(assignment)
  standardGeneric('iterations'))

#' @rdname accessors

setMethod('iterations',signature = 'Assignment',
          function(assignment){
            c(names(assignment@addIsoAssign),
              names(assignment@transAssign))
          })

#' @rdname accessors
#' @export

setGeneric('graph',function(assignment,
                            iteration,
                            type = c('filtered','all'))
  standardGeneric('graph'))

#' @rdname accessors

setMethod('graph',signature = 'Assignment',
          function(assignment, 
                   iteration,
                   type = c('filtered','all')){
            
            if (!iteration %in% iterations(assignment)) {
              iters <- assignment %>% 
                iterations() %>% 
                paste0('"',.,'"') %>% 
                paste(collapse = ', ')
              stop(paste0('Iteration not recognised. Argument `iteration` should be one of ',iters),
                   call. = FALSE)
            }
            
            type <- match.arg(type,
                              choices = c('filtered','all'))
            graph <- switch(type,
                            filtered = assignment@addIsoAssign[[iteration]]$filtered_graph,
                            all = assignment@addIsoAssign[[iteration]]$graph)
            
            return(graph)
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

#' assignedData
#' @rdname assignedData
#' @description Return data table used for assignments with feature assignments added to column names.
#' @param assignment S4 object of class Assignment
#' @return A tibble containing the original feature data with molecular formula assignments added to teh column names.
#' @examples 
#' \dontrun{
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' assignment <- assignMFs(feature_data,p)
#' 
#' assignedData(assignment)
#' }
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
#' @rdname summariseAssignments
#' @description Summarise features assigned to molecular formulas.
#' @param assignment S4 object of class Assignment
#' @return A tibble containing the feature assignments summarised by molecular formula.
#' @examples 
#' \dontrun{
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' assignment <- assignMFs(feature_data,p)
#' 
#' summariseAssignments(assignment)
#' }
#' @importFrom dplyr desc
#' @export

setGeneric('summariseAssignments',function(assignment)
  standardGeneric('summariseAssignments'))

#' @rdname summariseAssignments

setMethod('summariseAssignments',signature = 'Assignment',
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