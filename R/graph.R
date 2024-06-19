cond_node_key <- function(node) {
  
  cond_set <- node$conditional
  paste0(c("(", node$vertice, "| ",
           paste0(c(sapply(cond_set[-length(cond_set)],
                           function(cond_var) paste0(c(cond_var, ", "),
                                                     collapse = "")),
                    cond_set[length(cond_set)]),
                  collapse = ""),
           ")"), collapse = "")
}

CondNode <- function(vertice, conditional = c()) {
  
  node <- list(vertice = vertice,
               conditional=conditional)
  
  class(node) <- "CondNode"
  attr(node, "call") <- sys.call()
  
  node
  
}

print.CondNode <- function(x, ...) {
  #' @export
  cat(paste0(cond_node_key(x), "\n"))
}

print.Edge <- function(x, ...) {
  #' @export
  cat(paste0(cond_node_key(x[[1]]), " - ", cond_node_key(x[[2]]), "\n"))
}

equal_nodes <- function(n1, n2) {
  n1$vertice == n2$vertice & setequal(n1$conditional, n2$conditional)
}

NodeList <- function(nodes) {
  
  nodes <- nodes
  class(nodes) <- "NodeList"
  attr(nodes, "call") <- sys.call()
  
  nodes
}

print.NodeList <- function(x, all=FALSE) {
  #' @export
  if (all) {
    for (node in x){
      print(node)
    }
  } else {
    cat(paste0("NodeList contating ", length(x), " CondNode objects."))
  }
}

Edge <- function(n1, n2){
  edge <- list(n1, n2)
  class(edge) <- c("Edge", "NodeList")
  attr(edge, "call") <- sys.call()
  edge
}

edge_key <- function(edge) {
  paste0(cond_node_key(edge[[1]]), " - ", cond_node_key(edge[[2]]))
}

TreeGraph <- function(nodes, edges = NULL, pair_copulas = list()) {
  
  tree <- list(nodes = nodes, edges = edges,
               graph = igraph::make_empty_graph(n = length(nodes),
                                                directed = FALSE),
               pair_copulas = pair_copulas)
  
  class(tree) <- "TreeGraph"
  attr(tree, "call") <- sys.call()
  
  tree
  
}

TreeGraphList <- function(d) {
  trees <- list()
  trees <- append(
    trees, list(TreeGraph(NodeList(lapply(1:d, function(j) CondNode(j)))))
  )
  class(trees) <- "TreeGraphList"
  attr(trees, "call") <- sys.call()
  
  trees
}

all_initially_valid_edges <- function(nodes, which_include) {
  edges <- list()
  edge_mat <- NULL
  for (i in seq(1, length(nodes) - 1)) {
    for (j in seq(i + 1, length(nodes))) {
      if (valid_conditioning(nodes[[i]], nodes[[j]], which_include)) {
        edges <- append(edges, list(Edge(nodes[[i]], nodes[[j]])))
        if (is.null(edge_mat)) {
          edge_mat <- matrix(c(i, j), ncol = 1)
        } else {
          edge_mat <- cbind(edge_mat, c(i, j))
        }
      }
    }
  }
  list(edges = edges, edge_mat = edge_mat)
}

valid_conditioning <- function(n1, n2, which_include) {
  
  if (length(n1$conditional) != length(n2$conditional)){
    stop("Conditional sets of nodes must be of the same length!\n")
  }
  
  if (length(n1$vertice) != 1 | length(n2$vertice) != 1) {
    stop("A nodes vertice must contain exactly one element")
  }
  
  if (n1$vertice == n2$vertice) {
    FALSE
  } else if (length(intersect(n1$conditional, n2$conditional)) < 1 &
             length(n1$conditional) > 0) {
    FALSE
  } else if (n1$vertice %in% n2$conditional | n2$vertice %in% n1$conditional){
    FALSE
  } else if (!is.null(n1$conditional) &
             !(all(sort(n1$conditional) == sort(n2$conditional)))) {
    FALSE
  } else if (!(n1$vertice %in% which_include) | 
             !(n2$vertice %in% which_include)) {
    FALSE
  }
  else {
    TRUE
  }
  
}

add_edge_to_tree <- function(tree, edge) {
  tree$edges <- append(tree$edges, list(edge))
  tree
}

add_pair_copulas_to_tree <- function(tree, copula_pair) {
  tree$pair_copulas <- append(tree$pair_copulas, list(copula_pair))
  tree
}

edges_to_next_trees_nodes <- function(edges){
  
  new_nodes <- list()
  
  for (j in 1:length(edges)){
    
    if (length(condition_node(edges[[j]][[1]], edges[[j]][[2]])$conditional)
        == length(edges[[j]][[1]]$conditional) + 1) {
      
      new_nodes <- append(
        new_nodes, list(condition_node(edges[[j]][[1]],
                                       edges[[j]][[2]]))
      )
    }
    
    if (length(condition_node(edges[[j]][[2]], edges[[j]][[1]])$conditional)
        == length(edges[[j]][[2]]$conditional) + 1) {
      
      new_nodes <- append(
        new_nodes, list(condition_node(edges[[j]][[2]],
                                       edges[[j]][[1]]))
      )
    }
  }
  
  new_nodes
}

valid_edge <- function(tree, new_edge, new_edge_alias) {
  
  if (length(tree$edges) == 0) {
    TRUE
  } else {
    
    node1_edges <- sum(
      sapply(1:length(tree$edges),
             function(j) node_in_edge(new_edge[[1]], tree$edges[[j]]))
    )
    
    node2_edges <- sum(
      sapply(1:length(tree$edges),
             function(j) node_in_edge(new_edge[[2]], tree$edges[[j]]))
    )
    
    if (!(node1_edges > 0 & node2_edges > 0)) {
      TRUE
    } else {
      # Need to check if new edge makes a cycle.
      tmp_graph <- igraph::add.edges(tree$graph, new_edge_alias)
      
      igraph::girth(tmp_graph)$girth == 0
    }
  }
}

node_in_edge <- function(node, edge) {
  equal_nodes(node, edge[[1]]) | equal_nodes(node, edge[[2]])
}

condition_node <- function(conditioned, conditional){
  
  new_node <- list(
    vertice = conditioned$vertice,
    conditional = union(conditional$vertice,
                        union(conditioned$conditional,
                              conditional$conditional))
  )
  
  class(new_node) <- "CondNode"
  attr(new_node, "call") <- sys.call()
  
  new_node
}

equal_edges <- function(e1, e2) {
  ((equal_nodes(e1[[1]], e2[[1]]) & equal_nodes(e1[[2]], e2[[2]])) |
     (equal_nodes(e1[[1]], e2[[2]]) & equal_nodes(e1[[2]], e2[[1]])))
}
