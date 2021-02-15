library(DescTools)
library(igraph)
library(matrixStats)

build.assoc.graph = function(A, B, adjust = TRUE) {
  # Bonferroni adjusted confidence level
  if(adjust) alpha = 1 - 0.05 / 12
  else alpha = 1 - 0.05
  
  # flatted correlation matrices
  A_flat = t(sapply(A, unlist))
  B_flat = t(sapply(B, unlist))
  
  # average correlation per group
  A_means = colMeans(A_flat)
  B_means = colMeans(B_flat)
  
  # thresholds
  t = quantile(c(A_means, B_means), probs = .8, names = F)
  
  # Fisher confidence intervals
  A_fisher = sapply(A_means, function(x) CorCI(x, n = 12, conf.level = alpha))
  B_fisher = sapply(B_means, function(x) CorCI(x, n = 12, conf.level = alpha))
  
  # adjacency matrices
  n = dim(A_fisher)[2]
  A_matrix = sapply(1:n, function(i) int.test(A_fisher[2,i], A_fisher[3,i], t))
  B_matrix = sapply(1:n, function(i) int.test(B_fisher[2,i], B_fisher[3,i], t))
  A_matrix = as.integer(matrix(A_matrix, nrow = 116, byrow = TRUE))
  B_matrix = as.integer(matrix(B_matrix, nrow = 116, byrow = TRUE))
  
  return(list(A = A_matrix, B = B_matrix))
}

# correlation matrix for asd and td groups
cor.asd <- lapply(asd_sel, function(x) cor(x, method = "pearson"))
cor.td <- lapply(td_sel, function(x) cor(x, method = "pearson"))

# test for empty intersection of [-t,t] with CI(j,k)
# TRUE iff H_0 is rejected
int.test <- function(lower, upper, t) (t < lower) | (-t > upper)

# building adjacency matrix for ASD group
G_diff_bonferroni = build.assoc.graph(cor.asd, cor.td)
G_diff_raw = build.assoc.graph(cor.asd, cor.td, adjust = F)

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
t_boot = 0

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
  t_boot = t_boot + quantile(c(pooled_a, pooled_b), probs = 0.8, names = F)
  
  diff_graph_list[[i]] = abs(pooled_a - pooled_b)
}

# Normal 95% confidence interval for features correlations
flattened_features = matrix(unlist(diff_graph_list, recursive = TRUE), nrow = b)
t_boot = t_boot / b
means = colMeans(flattened_features)
deviations = colSds(flattened_features)
z = qnorm(1-0.05/12)
lower = matrix(means - z * deviations, nrow = 116, ncol = 116)
upper = matrix(means + z * deviations, nrow = 116, ncol = 116)
