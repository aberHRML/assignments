#' Extract graph attributes
#' @rdname graph
#' @description Extract node or edge attributes from a *tidygraph* `tbl_graph` object.
#' @param graph object of class tbl_graph
#' @examples 
#' a_graph <- tidygraph::tbl_graph(
#'   nodes = data.frame(
#'     name = c('a','b','c')
#'   ), 
#'   edges = data.frame(
#'     from = c(1,2),
#'     to = c(2,3),
#'     type = c(1,2)
#'   ))
#' 
#' ## Extract graph nodes
#' nodes(a_graph)
#' 
#' ## Extract graph edges
#' edges(a_graph)
#' @importFrom tibble as_tibble
#' @export

nodes <- function(graph){
  graph %>%
    vertex.attributes() %>%
    as_tibble()
}

#' @rdname graph
#' @importFrom igraph edge.attributes
#' @export

edges <- function(graph){
  graph %>%
    edge.attributes() %>%
    as_tibble()
}
