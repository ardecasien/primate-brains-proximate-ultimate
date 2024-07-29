my_plot_pcor <- function(M, row_names, edge_labels, thresh_nodes, thresh_edges, main) {
  G <- round(cov2pcor(M), digits = 3)
  N <- graph_from_adjacency_matrix(G,
                                   weighted = TRUE,
                                   mode = "undirected", diag = FALSE)
  
  layout_n <- graph_from_adjacency_matrix(abs(G),
                                          weighted = TRUE,
                                          mode = "undirected", diag = FALSE)
  
  my_palette <- colorRampPalette(c("red", "gray", "blue"))(n = 100)
  maxw <- max(E(N)$weight)
  
  E(layout_n)$color <- my_palette[1 + 99 * (E(N)$weight + maxw) / (2 * maxw)]
  E(layout_n)$edge.label <- edge_labels
  E(layout_n)$edge.label[edge_labels == ""] <- ""
  
  edge_width <- E(layout_n)$weight
  edge_width[edge_width < thresh_edges] <- 0
  
  plot(layout_n,
       rescale = TRUE,
       edge.width = 10 * edge_width,
       edge.label = E(layout_n)$edge.label,
       edge.label.family = "Helvetica",
       edge.label.color = "black",
       edge.label.cex = 0.7,
       vertex.label = row_names,
       vertex.label.family = "Helvetica",
       vertex.label.font = 2,
       vertex.label.cex = 0.8,
       vertex.color = "gray",
       vertex.label.color = "black",
       vertex.frame.color = "gray",
       main = main)
}

cov2pcor <- function(C) {
  C_inv <- solve(C)
  d <- sqrt(diag(C_inv))
  P <- -C_inv / outer(d, d)
  diag(P) <- 1
  return(P)
}

edge_exclusion_test <- function(C, n) {
  m <- cov2pcor(C)
  mt <- m[upper.tri(m)]
  edges <- -n * log(1 - mt^2)
  p_values <- pchisq(edges, df = 1, lower.tail = FALSE)
  significance_levels <- ifelse(p_values < 0.05, "*","")
  
  # Create a matrix with empty strings
  result <- matrix("", nrow = nrow(m), ncol = ncol(m))
  
  # Assign significance levels to upper triangular part
  result[upper.tri(result)] <- significance_levels
  
  # Return the result matrix and p-values
  return(list(matrix = result, p_values = p_values))
}
