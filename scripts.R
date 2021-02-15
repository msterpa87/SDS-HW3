library(DescTools)
library(igraph)
library(matrixStats)

build_graph = function(A, B, probs = 0.8, adjust = TRUE) {
  # Bonferroni adjusted confidence level
  if(adjust) alpha = 1 - 0.05 / 12
  else alpha = 1 - 0.05
  
  # correlation matrices
  A = lapply(A, function(x) cor(x, method = "pearson"))
  B = lapply(B, function(x) cor(x, method = "pearson"))
  
  # flatted correlation matrices
  A_flat = t(sapply(A, unlist))
  B_flat = t(sapply(B, unlist))
  
  # thresholds
  t_a = quantile(A_flat, probs = probs, names = F)
  t_b = quantile(B_flat, probs = probs, names = F)
  
  # average correlation per group
  A_means = colMeans(A_flat)
  B_means = colMeans(B_flat)
  
  # Fisher confidence intervals
  A_fisher = sapply(A_means, function(x) CorCI(x, n = 145, conf.level = alpha))
  B_fisher = sapply(B_means, function(x) CorCI(x, n = 145, conf.level = alpha))
  
  # adjacency matrices
  n = dim(A_fisher)[2]
  A_matrix = sapply(1:n, function(i) int.test(A_fisher[2,i], A_fisher[3,i], t_a))
  B_matrix = sapply(1:n, function(i) int.test(B_fisher[2,i], B_fisher[3,i], t_b))
  A_matrix = matrix(as.integer(A_matrix), nrow = 116, byrow = TRUE)
  B_matrix = matrix(as.integer(B_matrix), nrow = 116, byrow = TRUE)
  
  return(list(A = A_matrix, B = B_matrix))
}

plot_graph = function(M, title) {
  G <- simplify(graph_from_adjacency_matrix(M), remove.loops = T)
  V(G)$size <- 1
  V(G)$label <- ""
  E(G)$arrow.mode <- 0
  plot(G, main = title)
}


# test for empty intersection of [-t,t] with CI(j,k)
# TRUE iff H_0 is rejected
int.test <- function(lower, upper, t) (t < lower) | (-t > upper)

# building adjacency matrix for ASD group
G_bonferroni = build_graph(asd_sel, td_sel)
G_raw = build_graph(asd_sel, td_sel, adjust = F)

# node degree histogram in difference graph
par(mfrow = c(1,2))
hist(rowSums(G_diff_bonferroni), col = "orchid", border = "white",
     main = "Bonferroni adjusted", xlab = "Node degree")
hist(rowSums(G_diff_raw), col = "pink", border = "white",
     main = "Non adjusted", xlab = "Node degree")

# Plotting graphs
par(mfrow = c(2,2))
plot_graph(G_bonferroni$A, "Group A (Bonferroni)")
plot_graph(G_bonferroni$B, "Group B (Bonferroni)")
plot_graph(G_raw$A, "Group A (non-adjusted)")
plot_graph(G_raw$B, "Group B (non-adjusted)")

# histograms of co-activation with sliding threshold
probs = seq(0.5, 1, 0.05)
n = length(probs)
n_coact_bonf = rep(NA, n)
n_coact_raw = rep(NA, n)

for(i in 1:n) {
  G_bonf = build_graph(asd_sel, td_sel, probs = probs[i], adjust = T)
  G_raw = build_graph(asd_sel, td_sel, probs = probs[i], adjust = F)
  
  # number of paired co-activation between groups
  n_coact_bonf[i] = sum((G_bonf$A == 1) & (G_bonf$B == 1)) - 116
  n_coact_raw[i] = sum((G_raw$A == 1) & (G_raw$B == 1)) - 116
}

par(mfrow = c(1,1))
plot(x = probs, y = n_coact_bonf, type = "b", col = "red", pch = 16, xlab = "Threshold",
     ylab = "Number of Paired Co-Activations", ylim = range(n_coact_raw),
     main = "Co-activation at different thresholds")
lines(x = probs, y = n_coact_raw, type = "b", col = "blue", pch = 16, add = TRUE)
segments(.5, n_coact[7], 1, n_coact[7], lwd = 2, lty = 2, col = "green")
text(.8, n_coact[7] - 150, bquote(q[80]))
legend("topright", legend=c("Bonferroni", "Non-Adjusted"),
       col=c("red", "blue"), lty=1, lwd=3)
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
