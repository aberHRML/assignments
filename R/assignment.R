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
#' @param component component number to extract
#' @param feature feature information to extract
#' @examples 
#' mf_assignments <- new('Assignment',
#'                   data = feature_data)
#' 
#' ## Return feature data
#' featureData(mf_assignments)
#' 
#' ## Return correlations
#' correlations(mf_assignments)
#'
#' ## Return relationships
#' relationships(mf_assignments)
#' 
#' ## Return the available iterations
#' iterations(mf_assignments)
#'
#' ## Return a selected graph
#' \dontrun{
#' graph(mf_assignments,'A&I1')
#' }
#'
#' ## Return a component information for a selected graph
#' \dontrun{
#' components(mf_assignments,'A&I1')
#' }
#' 
#' ## Return a component information for a selected feature
#' \dontrun{
#' featureComponents(mf_assignments,'n191.01962')
#' }
#'
#'  ## Extract a component graph
#' \dontrun{
#' component(mf_assignments,1,'A&I1')
#' }
#' 
#' ## Return assignments
#' assignments(mf_assignments)
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
            assignment@correlations
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
                            type = c('selected','all'))
  standardGeneric('graph'))

#' @rdname accessors

setMethod('graph',signature = 'Assignment',
          function(assignment, 
                   iteration,
                   type = c('selected','all')){
            
            if (!iteration %in% iterations(assignment)) {
              iters <- assignment %>% 
                iterations() %>% 
                paste0('"',.,'"') %>% 
                paste(collapse = ', ')
              stop(paste0('Iteration not recognised. Argument `iteration` should be one of ',iters),
                   call. = FALSE)
            }
            
            type <- match.arg(type,
                              choices = c('selected','all'))
            
            assignment_iteration <- switch(
              str_remove(iteration,'[1-9]'),
              `A&I` = assignment@addIsoAssign,
              `T` = assignment@transAssign)
            
            graph <- switch(type,
                            selected = assignment_iteration[[iteration]]$filtered_graph,
                            all = assignment_iteration[[iteration]]$graph)
            
            return(graph)
          })

#' @rdname accessors
#' @export

setGeneric('components',function(assignment,
                                 iteration,
                                 type = c('selected','all'))
  standardGeneric('components'))

#' @rdname accessors

setMethod('components',signature = 'Assignment',
          function(assignment, 
                   iteration,
                   type = c('selected','all')){
            
            selected_graph <- graph(assignment,iteration,type) %>% 
              nodes() %>% 
              select(Component:`Component Plausibility`) %>% 
              distinct()
            
            return(selected_graph)
          })

#' @rdname accessors
#' @export

setGeneric('featureComponents',function(assignment,
                                        feature,
                                        type = c('selected','all'))
  standardGeneric('featureComponents'))

#' @rdname accessors

setMethod('featureComponents',signature = 'Assignment',
          function(assignment, 
                   feature,
                   type = c('selected','all')){
            
            available_iterations <- iterations(assignment)
            
            available_iterations %>% 
              map(graph,assignment = assignment,type = type) %>% 
              set_names(available_iterations) %>% 
              map_dfr(nodes,.id = 'Iteration') %>% 
              filter(Feature == feature)
          })

#' @rdname accessors
#' @export

setGeneric('component',function(assignment,
                                component,
                                iteration,
                                type = c('selected','all'))
  standardGeneric('component'))

#' @rdname accessors

setMethod('component',signature = 'Assignment',
          function(assignment, 
                   component,
                   iteration,
                   type = c('selected','all')){
            
            iteration_components <- components(assignment,
                                               iteration,
                                               type)
            
            if (!component %in% iteration_components$Component){
              stop(paste0('Component ',component, ' not found in iteration ',iteration,'.'))
            }
            
            graph(assignment,iteration,type) %>% 
              filter(Component == component)
            
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
              featureData()
            
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
#' mf_assignments <- assignMFs(feature_data,p)
#' 
#' summariseAssignments(mf_assignments)
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

#' Create an Assignment S4 class object
#' @rdname assignment
#' @description Constructor methods for creating an object of S4 class `Assignment`.
#' @param feature_data a tibble or an object of S4 class `AnalysisData` or `Analysis` containing the feature intensity matrix of m/z for which to assign molecular formulas. See details.
#' @param parameters an S4 object of class `AssignmentParamters` containing the parameters for molecular formula assignment
#' @return An object of S4 class `Assignment`.
#' @examples 
#' mf_assignments <- assignment(feature_data,assignmentParameters('FIE-HRMS'))
#' @export

setGeneric('assignment',function(feature_data,parameters,...)
  standardGeneric('assignment'))

#' @rdname assignment

setMethod('assignment',signature = c('tbl_df','AssignmentParameters'),
          function(feature_data,parameters){
            new('Assignment',
                parameters,
                data = feature_data)
          })


#' @rdname assignment

setMethod('assignment',signature = c('AnalysisData','AssignmentParameters'),
          function(feature_data,parameters){
            new('Assignment',
                parameters,
                data = feature_data%>% 
                  dat())
          })

#' @rdname assignment

setMethod('assignment',signature = c('Analysis','AssignmentParameters'),
          function(feature_data,parameters,type = c('pre-treated','raw')){
            
            type <- match.arg(type,
                              choices = c('pre-treated','raw'))
            
            if (type == 'raw') feature_data <- raw(feature_data)
            if (type == 'pre-treated') feature_data <- preTreated(feature_data)
            
           assignment(feature_data,
                      parameters)
          })