library(DescTools)
library(igraph)
library(matrixStats)

# standardizing data
asd_sel = scale(asd_sel)
td_sel = scale(td_sel)

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

# threshold (80th percentile) of the combined correlations of both groups
cor.all <- c(unlist(cor.asd.pooled), unlist(cor.td.pooled))
threshold <- as.numeric(quantile(cor.all, probs = 0.8))

# test for empty intersection of [-t,t] with CI(j,k)
# TRUE iff H_0 is rejected
int.test <- function(lower, upper, t) (t < lower) | (-t > upper)

build.assoc.graph <- function(cor.matrices.a, cor.matrices.b, threshold, D = 116, m = 12,
                              bonferroni = TRUE) {
  # D : features, n: time series length, m: patients
  adj.matrix <- matrix(nrow = D, ncol = D)
  
  # bonferroni adjusted confidence level
  if(bonferroni) alpha = 1 - 0.05 / m
  else alpha = 1 - 0.05
  
  # building adjacency matrix
  for(j in 1:D) {
    for(k in 1:D) {
      
      # vector of pearson correlation of features (j,k) for each patient
      p1 <- as.numeric(unlist(lapply(cor.matrices.a, function(x) x[j, k])))
      p2 <- as.numeric(unlist(lapply(cor.matrices.b, function(x) x[j, k])))
      p <- abs(p1 - p2)

      # Fisher CI for each patient correlation coefficient
      conf.int <- lapply(p, function(x) CorCI(x, n = m, conf.level = alpha))
    
      # conf.int[i] = (rho, lower, upper)
      # bool vector of CI intersection test with threshold interval
      int.test.results <- unlist(lapply(conf.int, function(x) int.test(x[2], x[3], threshold)))
      
      # edge (j,k) iff at least one H_0 is rejected
      adj.matrix[j,k] <- as.numeric(any(int.test.results, na.rm = T))
    }
  }
  return(adj.matrix)
}

# building adjacency matrix for ASD group
G_diff_bonferroni = build.assoc.graph(cor.asd, cor.td, threshold)
G_diff_raw = build.assoc.graph(cor.asd, cor.td, threshold, bonferroni = F)

# node degree histogram in difference graph
par(mfrow = c(1,2))
hist(rowSums(G_diff_bonferroni), col = "orchid", border = "white",
     main = "Bonferroni adjusted", xlab = "Node degree")
hist(rowSums(G_diff_raw), col = "pink", border = "white",
     main = "Non adjusted", xlab = "Node degree")

# simple plotting functions
graph.asd <- graph_from_adjacency_matrix(asd.adj.matrix)
V(graph.asd)$size <- 1
V(graph.asd)$label <- ""
E(graph.asd)$arrow.mode <- 0
plot(graph.asd)

# number of paired co-activation between groups at different thresholds
A = cor.asd.pooled
B = cor.td.pooled
t_vec = quantile(c(A, B), seq(0.5, 1, 0.05), names = F)
y = unlist(lapply(t_vec, function(t) sum((A >= t) & (B >= t))))
plot(x = t_vec, y, col = "orchid", pch = 16, xlab = "Threshold",
     ylab = "Number of Paired Co-Activations",
     main = "Co-Activations between groups at different thresholds")
segments(0, y[6], 1, y[6], lwd = 2, lty = 2, col = "green")
text(0.15, 2800, "80th percentile")
grid()

# Bootstrap
b = 10
diff_graph_list = list()

for(i in 1:b) {
  # resample of each group
  sample_a = sample(asd_sel, 12, replace = T)
  sample_b = sample(td_sel, 12, replace = T)
  
  # correlation matrix for asd and td groups
  grouped_sample_a = lapply(sample_a, function(x) cor(x, method = "pearson"))
  grouped_sample_b = lapply(sample_b, function(x) cor(x, method = "pearson"))
  
  # pooling taking the median
  pooled_a = pooling.cor(grouped_sample_a)
  pooled_b = pooling.cor(grouped_sample_b)
  
  diff_graph_list[[i]] = abs(pooled_a - pooled_b)
}

means = colMeans(flattened_features)
deviations = colSds(flattened_features)
lower = matrix(means - z * deviations, nrow = 116, ncol = 116)
upper = matrix(means + z * deviations, nrow = 116, ncol = 116)
