#' assignments-Assignment
#' @rdname assignments
#' @description Get table of assigned features from an Assignment
#' @param assignment S4 object of class Assignment
#' @export

setMethod('assignments',signature = 'Assignment',
          function(assignment){
            assignment@assignments
})

#' nodes
#' @description extract node table from tbl_graph object.
#' @param graph object of class tbl_graph
#' @export

nodes <- function(graph){
           graph %>%
              vertex.attributes() %>%
              as_tibble()
}

#' edges
#' @description extract edge table from tbl_graph object.
#' @param graph object of class tbl_graph
#' @importFrom igraph edge.attributes
#' @export

edges <- function(graph){
  graph %>%
    edge.attributes() %>%
    as_tibble()
}

#' assignmentData
#' @rdname assignmentData
#' @description Return data table used for assignments.
#' @param assignment S4 object of class Assignment
#' @export

setMethod('assignmentData', signature = 'Assignment',
          function(assignment){
            assignment@data
})

#' assignedData
#' @rdname assignedData
#' @description Return data table used for assignments with feature assignments added to column names.
#' @param assignment S4 object of class Assignment
#' @export

setMethod('assignedData', signature = 'Assignment',
          function(assignment){
            
            d <- x %>%
              assignmentData()
            
            assignedFeats <- x %>%
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