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
