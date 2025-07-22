# Define functions for doing the random acquisition simulation
multisample <- function(n, x, fct, simplify2vec = FALSE) {
  y <- mapply(sample, split(x, fct), n, SIMPLIFY = FALSE)
  return(if(simplify2vec) unname(do.call(c, y)) else y)
}

balanced_ranet <- function(vocab_size, G, POS) {
  env_size <- igraph::vcount(G)
  assertthat::are_equal(env_size, length(POS))
  assertthat::are_equal(nlevels(POS), length(vocab_size))
  ix <- multiSample(vocab_size, seq_len(env_size), POS, simplify2vec = TRUE)
  return(igraph::induced_subgraph(G, ix))
}

network_stats <- function(g) {
    return(c(
        indegree = median(igraph::degree(g, mode = "in")),
        clustcoef = igraph::transitivity(g, type = "global"),
        meandist = igraph::mean_distance(g))
    )
}

# If we run this function, each random graph will be discarded after computing
# statistics. This eliminates the memory consumption issue.
balanced_ranet_stats <- function(vocab_size, g, pos) {
    return(network_stats(balanced_RAN_network(vocab_size, g, pos)))
}

