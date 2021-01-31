# correlation matrix for asd and td groups
cor.asd <- lapply(asd_sel, cor)
cor.td <- lapply(td_sel, cor)

# combined correlations from both groups
cor.all <- c(unlist(cor.asd), unlist(cor.td))

# threshold (80th percentile)
t <- quantile(cor.all, probs = 0.8)

# test for empty intersection of [-t,t] with CI(j,k)
# TRUE iff H_0 is rejected
int.test <- function(lower, upper, t) (-t < lower | -t > upper) & (t < lower | t > upper)

build.assoc.graph <- function(cor.matrices, D = 116, n = 145, m = 12) {
  # D : features, n: time series length, m: patients
  adj.matrix <- matrix(nrow = D, ncol = D)
  
  # bonferroni adjusted confidence level
  alpha <- 1 - 0.05 / m
  
  # building adjacency matrix
  for(j in 1:D) {
    for(k in 1:D) {
      
      # vector of pearson correlation of features (j,k) for each patient
      p <- as.numeric(unlist(lapply(cor.matrices, function(x) x[j, k])))
      
      # Fisher CI for each patient corrleation coefficient
      conf.int <- lapply(p, function(x) CorCI(x, n = n, conf.level = alpha))
    
      # conf.int[i] = (rho, lower, upper)
      # bool vector of CI intersection test with threshold interval
      int.test.results <- unlist(lapply(conf.int, function(x) int.test(x[2], x[3], t)))
      
      # edge (j,k) iff at least one H_0 is rejected
      adj.matrix[j,k] <- as.numeric(any(int.test.results))
    }
  }
  return(adj.matrix)
}

# building adjacency matrix for ASD group
asd.adj.matrix <- build.assoc.graph(cor.asd)

# simple plotting functions
graph.asd <- graph_from_adjacency_matrix(asd.adj.matrix)
V(graph.asd)$size <- 1
V(graph.asd)$label <- ""
E(graph.asd)$arrow.mode <- 0
plot(graph.asd)
