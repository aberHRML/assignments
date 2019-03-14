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
#' @export

nodes <- function(graph){
           graph %>%
              vertex.attributes() %>%
              as_tibble()
}

#' edges
#' @importFrom igraph edge.attributes
#' @export

edges <- function(graph){
  graph %>%
    edge.attributes() %>%
    as_tibble()
}