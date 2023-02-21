#' Assignment
#' @rdname Assignment-class
#' @description An S4 class to store molecular formula assignment results.
#' @slot log a list containing assignment logs
#' @slot data a tibble containing the *m/z* peak intensity matrix
#' @slot correlations a tibble containing the correlations analysis results
#' @slot relationships a tibble containing the calculated mathematical relationships
#' @slot addIsoAssign a list containing the results of the adduct and isotope assignment iterations
#' @slot transAssign a list containing the results of the transformation assignment iterationst
#' @slot assignments a tibble containing the assigned molecular formulas

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
#' @description Access methods for `Assignment` S4 class
#' @param assignment S4 object of class Assignment
#' @param iteration the assignment iteration
#' @param type the graph type to return. `filtered` returns the assignment graph after component selection. `all` returns all assignment components.
#' @param component component number to extract
#' @param feature feature information to extract
#' @details 
#' * `featureData` - Return the initially specifed *m/z* feature data.
#' * `correlations` - Return the correlation analysis results.
#' * `relationships` - Return the calculated relationships.
#' * `iterations` - Return the assignment iteration performed.
#' * `graph` - Return a selected graph.
#' * `components` - Return the component information for an assignment iteration.
#' * `featureComponents` - Return the component information for a selected feature.
#' * `component` - Extract a component graph.
#' * `assignments` - Return the molecular formulas assigned to the *m/z* features.
#' * `assignedData` - Return the *m/z* peak intensity matrix with the molecular formula assignments included in the column names.
#' * `summariseAssignments` - Return a tibble of the assignments summarised by molecular formula.
#' @return A tibble or `tbl_graph` containing assignment results depending on the method used. 
#' @examples 
#' plan(future::sequential)
#' p <- assignmentParameters('FIE-HRMS')
#'
#' mf_assignments <- assignMFs(feature_data,p)
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
#' graph(mf_assignments,'A&I1')
#'
#' ## Return a component information for a selected graph
#' components(mf_assignments,'A&I1')
#' 
#' ## Return a component information for a selected feature
#' featureComponents(mf_assignments,'n191.01962')
#'
#'  ## Extract a component graph
#' component(mf_assignments,1,'A&I1')
#' 
#' ## Return assignments
#' assignments(mf_assignments)
#' 
#' ## Return an m/z intensity matrix with the assignments included
#' ## in the column names
#' assignedData(mf_assignments)
#' 
#' ## Return the assignments summarised by molecular formula
#' summariseAssignments(mf_assignments)
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

#' @rdname accessors
#' @export

setGeneric('assignedData',function(assignment)
  standardGeneric('assignedData'))

#' @rdname accessors

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

#' @rdname accessors
#' @importFrom dplyr desc
#' @export

setGeneric('summariseAssignments',function(assignment)
  standardGeneric('summariseAssignments'))

#' @rdname accessors

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
#' @param type type `pre-treated` or `raw` data on which to perform assignment when argument `feature_data` is of class `Analysis`
#' @param ... arguments to pass to the relevant method
#' @return An object of S4 class `Assignment`.
#' @examples 
#' mf_assignments <- assignment(feature_data,assignmentParameters('FIE-HRMS'))
#' mf_assignments
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