# correlation matrix for asd and td groups
cor.asd <- lapply(asd_sel, function(x) cor(x, method = "pearson"))
cor.td <- lapply(td_sel, function(x) cor(x, method = "pearson"))

pooling.cor <- function(cor.matrices, D = 116) {
  M <- matrix(nrow = D, ncol = D)
  for(i in 1:D) {
    for(j in 1:D) {
      M[i,j] <- median(unlist(lapply(1:12, function(k) cor.matrices[[k]][i,j])))
    }
  }
  return(M)
}

cor.asd.pooled <- pooling.cor(cor.asd)
cor.td.pooled <- pooling.cor(cor.td)

# combined correlations from both groups
cor.all <- c(unlist(cor.asd.pooled), unlist(cor.td.pooled))

# threshold (80th percentile)
threshold <- as.numeric(quantile(cor.all, probs = 0.8))

# test for empty intersection of [-t,t] with CI(j,k)
# TRUE iff H_0 is rejected
int.test <- function(lower, upper, t) (t < lower) | (-t > upper)

build.assoc.graph <- function(cor.matrices.a, cor.matrices.b, threshold, D = 116, n = 145, m = 12) {
  # D : features, n: time series length, m: patients
  adj.matrix <- matrix(nrow = D, ncol = D)
  
  # bonferroni adjusted confidence level
  alpha <- 1 - 0.05 / m
  
  # building adjacency matrix
  for(j in 1:D) {
    for(k in 1:D) {
      
      # vector of pearson correlation of features (j,k) for each patient
      p1 <- as.numeric(unlist(lapply(cor.matrices.a, function(x) x[j, k])))
      p2 <- as.numeric(unlist(lapply(cor.matrices.b, function(x) x[j, k])))
      p <- abs(p1 - p2)
      print(p)
      
      # Fisher CI for each patient corrleation coefficient
      conf.int <- lapply(p, function(x) CorCI(x, n = n, conf.level = alpha))
    
      # conf.int[i] = (rho, lower, upper)
      # bool vector of CI intersection test with threshold interval
      int.test.results <- unlist(lapply(conf.int, function(x) int.test(x[2], x[3], threshold)))
      
      # edge (j,k) iff at least one H_0 is rejected
      adj.matrix[j,k] <- as.numeric(any(int.test.results))
    }
  }
  return(adj.matrix)
}

# building adjacency matrix for ASD group
asd.adj.matrix <- build.assoc.graph(cor.asd, cor.td, threshold)

# simple plotting functions
graph.asd <- graph_from_adjacency_matrix(asd.adj.matrix)
V(graph.asd)$size <- 1
V(graph.asd)$label <- ""
E(graph.asd)$arrow.mode <- 0
plot(graph.asd)
