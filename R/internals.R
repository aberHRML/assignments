
networkComponents <- function(nodes,edges){
  graph <- as_tbl_graph(edges,directed = F)
}

degree <- function(n_nodes,n_edges){
  2 * (n_edges / n_nodes)
}

plausibility <- function(AIS,degree,weight){
  AIS + weight + degree/10
}
