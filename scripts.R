library(DescTools)
library(igraph)
library(matrixStats)

load("C:/Users/Stefania/Desktop/Sapienza/I SEMESTRE_20-21/STATISTICAL METHODS IN DATA SCIENCE AND LABORATORY I/HW3/hw3_data.RData")

build_graph = function(A_flat, B_flat, probs = 0.8, adjust = TRUE) {
  # Bonferroni adjusted confidence level
  if(adjust) alpha = 1 - 0.05 / 12
  else alpha = 1 - 0.05
  
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
  V(G)$size <- .5
  V(G)$label <- ""
  E(G)$width = 0.01
  E(G)$arrow.mode <- 0
  plot(G, main = title, vertex.color="indianred3", vertex.size=5, vertex.frame.color="ivory1", 
       vertex.label.cex=0.8, vertex.label.dist=5, edge.curved=0.2)
}

plot_comparison_graphs = function(G_bonferroni, G_raw) {
  # Plotting graphs
  par(mfrow = c(2,2), mai=c(.1,.1,.2,.1))
  plot_graph(G_bonferroni$A, "Group A (Bonferroni)")
  plot_graph(G_bonferroni$B, "Group B (Bonferroni)")
  plot_graph(G_raw$A, "Group A (non-adjusted)")
  plot_graph(G_raw$B, "Group B (non-adjusted)")
}

# test for empty intersection of [-t,t] with CI(j,k)
# TRUE iff H_0 is rejected
int.test <- function(lower, upper, t) (t < lower) | (-t > upper)

# correlation matrices
A = lapply(asd_sel, function(x) cor(x, method = "pearson"))
B = lapply(td_sel, function(x) cor(x, method = "pearson"))

# flatted correlation matrices
A_flat = t(sapply(A, unlist))
B_flat = t(sapply(B, unlist))

# building adjacency matrix for ASD group
G_bonferroni = build_graph(A_flat, B_flat)
G_raw = build_graph(A_flat, B_flat, adjust = F)

plot_comparison_graphs(G_bonferroni, G_raw)

# histograms of co-activation with sliding threshold
probs = seq(0.5, 1, 0.05)
n = length(probs)
n_coact_bonf = rep(NA, n)
n_coact_raw = rep(NA, n)

for(i in 1:n) {
  G_bonf = build_graph(A_flat, B_flat, probs = probs[i], adjust = T)
  G_raw = build_graph(A_flat, B_flat, probs = probs[i], adjust = F)
  
  # number of paired co-activation between groups
  n_coact_bonf[i] = sum((G_bonf$A == 1) & (G_bonf$B == 1)) - 116
  n_coact_raw[i] = sum((G_raw$A == 1) & (G_raw$B == 1)) - 116
}

par(mfrow = c(1,1))
plot(x = probs, y = n_coact_bonf, type = "b", col = "seagreen3", pch = 16, xlab = "Threshold",
     ylab = "Number of Paired Co-Activations", ylim = range(n_coact_raw), lwd = 2,
     main = "Co-activation at different thresholds")
lines(x = probs, y = n_coact_raw, type = "b", col = "indianred3", pch = 16, lwd = 2,add = TRUE)
segments(.5, n_coact_bonf[7], 1, n_coact_bonf[7], lwd = 2, lty = 2, col = "lightcyan4")
text(.8, n_coact_bonf[7] - 150, bquote(q[80]),lwd = 2)
legend("topright", legend=c("Bonferroni", "Non-Adjusted"),
       col=c("seagreen3", "indianred3"), lty=1, lwd=3)
grid()

# Bootstrap
b = 10
A_flat_boot = matrix(0, nrow = 12, ncol = 116^2)
B_flat_boot = matrix(0, nrow = 12, ncol = 116^2)

for(i in 1:b) {
  # resample of each group
  A = sample(asd_sel, 12, replace = T)
  B = sample(td_sel, 12, replace = T)
  
  # correlation matrix for asd and td groups
  A = lapply(sample_a, function(x) cor(x, method = "pearson"))
  B = lapply(sample_b, function(x) cor(x, method = "pearson"))
  
  # flatted correlation matrices
  A_flat = t(sapply(A, unlist))
  B_flat = t(sapply(B, unlist))
  
  # average correlation per group
  A_flat_boot = A_flat_boot + A_flat
  B_flat_boot = B_flat_boot + B_flat
}

A_flat_boot = A_flat_boot / b
B_flat_boot = B_flat_boot / b

# building adjacency matrix for ASD group
G_bonferroni_boot = build_graph(A_flat_boot, B_flat_boot)
G_raw_boot = build_graph(A_flat_boot, B_flat_boot, adjust = F)

plot_comparison_graphs(G_bonferroni_boot, G_raw_boot)

# Difference graph
A_means = colMeans(A_flat_boot)
B_means = colMeans(B_flat_boot)
G_delta = abs(A_means - B_means)

# Plot of differences at different thresholds
t_range = quantile(G_delta, probs = c(0.85, 0.95))
t_diff = seq(t_range[1], t_range[2], length.out = 4)
adj_matrix = function(t) as.integer(G_delta >= t)
delta_matrices = t(sapply(t_diff, adj_matrix))
to_matrix = function(M_flat) matrix(M_flat, nrow=116, byrow=T)
titles = unlist(lapply(t_diff, function(t) paste(c("t =", round(t,2)), collapse = " ")))

plot_degree = function(M, title, color) {
  G = graph_from_adjacency_matrix(M)
  deg = degree(G, mode="all")
  deg.dist = degree.distribution(G, cumulative=T, mode="all")
  plot(x=0:max(deg), y=1-deg.dist, pch=19, type="l", col=color,
       xlab="Degree", ylab="Cumulative frequency", main=title,
       lwd=2)
  grid()
}

plot_degree_distributions = function(flat_matrices, titles) {
  par(mfrow = c(2,2), mai=c(.3,.3,.3,.3))
  colors = c("pink", "orange", "blue", "orchid")
  invisible(sapply(1:4, function(i) plot_degree(to_matrix(flat_matrices[i,]),
                                                titles[i], colors[i])))
}

plot_degree_distributions(delta_matrices, titles)

