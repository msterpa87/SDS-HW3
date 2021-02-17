library(DescTools)
library(igraph)
library(matrixStats)

load("C:/Users/Stefania/Desktop/Sapienza/I SEMESTRE_20-21/STATISTICAL METHODS IN DATA SCIENCE AND LABORATORY I/HW3/hw3_data.RData")

average_cor = function(M) FisherZInv(colMeans(FisherZ(M)))
to_matrix = function(M_flat) matrix(as.integer(M_flat), nrow=116, byrow=T)

build_graph = function(A_flat, B_flat, probs = 0.8, adjust = TRUE) {
  # Bonferroni adjusted confidence level
  d = 116
  m = d*(d-1)/2
  if(adjust) alpha = 1 - 0.05 / m
  else alpha = 1 - 0.05
  
  # average correlation per group
  A_means = average_cor(A_flat)
  B_means = average_cor(B_flat)
  
  # thresholds
  t_a = quantile(A_means, probs = probs, names = F, na.rm = T)
  t_b = quantile(B_means, probs = probs, names = F, na.rm = T)
  
  # Fisher confidence intervals
  A_fisher = sapply(A_means, function(x) CorCI(x, n = 145, conf.level = alpha))
  B_fisher = sapply(B_means, function(x) CorCI(x, n = 145, conf.level = alpha))
  
  # adjacency matrices
  n = dim(A_fisher)[2]
  A_matrix = sapply(1:n, function(i) int.test(A_fisher[2,i], A_fisher[3,i], t_a))
  B_matrix = sapply(1:n, function(i) int.test(B_fisher[2,i], B_fisher[3,i], t_b))
  A_matrix = to_matrix(A_matrix)
  B_matrix = to_matrix(B_matrix)
  
  return(list(A = A_matrix, B = B_matrix))
}

set_graph_params = function(G) {
  V(G)$size <- .5
  V(G)$label <- ""
  E(G)$width = 0.01
  E(G)$arrow.mode <- 0
  return(G)
}

plot_graph = function(M, title, community=F) {
  G = simplify(graph_from_adjacency_matrix(M, mode="undirected"), remove.loops = T)
  G = set_graph_params(G)
  params = list(G, main = title, vertex.color="indianred3", vertex.size=5, vertex.frame.color="ivory1", 
                vertex.label.cex=0.8, vertex.label.dist=5, edge.curved=0.2)
  
  # Plot with community detection
  if(community) {
    cfg = cluster_fast_greedy(G)
    params = append(list(cfg), params)
  }
  do.call(plot, params)
}

plot_comparison_graphs = function(G_bonferroni, G_raw, community=F) {
  # Plotting graphs
  par(mfrow = c(2,2), mai=c(.1,.1,.2,.1))
  graph_list = list(G_bonferroni$A, G_bonferroni$B, G_raw$A, G_raw$B)
  title_list = c("Group A (Bonferroni)", "Group B (Bonferroni)",
                 "Group A (non-adjusted)", "Group B (non-adjusted)")
  invisible(sapply(1:4, function(i) plot_graph(graph_list[[i]], title_list[[i]],
                                               community=community)))
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
plot_comparison_graphs(G_bonferroni, G_raw, community = T)

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

par(mfrow = c(1,1), mai=c(.5,.5,.5,.5))
plot(x = probs, y = n_coact_bonf, type = "b", col = "seagreen3", pch = 16, xlab = "Threshold",
     ylab = "Number of Paired Co-Activations", ylim = range(n_coact_raw), lwd = 2,
     main = "Co-activation at different thresholds")
lines(x = probs, y = n_coact_raw, type = "b", col = "indianred3", pch = 16, lwd = 2)
segments(.5, n_coact_bonf[7], 1, n_coact_bonf[7], lwd = 2, lty = 2, col = "lightcyan4")
text(.8, n_coact_bonf[7] - 150, bquote(q[80]),lwd = 2)
legend("topright", legend=c("Bonferroni", "Non-Adjusted"),
       col=c("seagreen3", "indianred3"), lty=1, lwd=3)
grid()

#######################################
#               BOOTSTRAP             #
#######################################

b = 100

# point estimate correlation
mu_cor_a = average_cor(A_flat)
mu_cor_b = average_cor(B_flat)

t_a = quantile(mu_cor_a, probs = .8, na.rm = T)
t_b = quantile(mu_cor_b, probs = .8, na.rm = T)

A_boot = matrix(NA, nrow = b, ncol = 116^2)
B_boot = matrix(NA, nrow = b, ncol = 116^2)

n = 12

for(i in 1:b) {
  idx_a = sample(1:n, n, replace = T)
  idx_b = sample(1:n, n, replace = T)
  
  A_boot[i,] = average_cor(A_flat[idx_a,])
  B_boot[i,] = average_cor(B_flat[idx_b,])
}

d = 116
m = d*(d-1)/2
z = qnorm(1-0.05/m)

# Bootstrap standard error
se_boot_a = apply(A_boot, MARGIN = 2, sd)
se_boot_b = apply(B_boot, MARGIN = 2, sd)

# Normal confidence intervals
lower_a = mu_cor_a - z * se_boot_a
upper_a = mu_cor_a + z * se_boot_a

lower_b = mu_cor_b - z * se_boot_b
upper_b = mu_cor_b + z * se_boot_b

# Adjacency matrices
M_a = to_matrix(sapply(1:116^2, function(i) int.test(lower_a[i], upper_a[i], t_a)))
M_b = to_matrix(sapply(1:116^2, function(i) int.test(lower_b[i], upper_b[i], t_b)))

par(mfrow=c(2,2), mai=c(.1,.1,.2,.1))
plot_graph(M_a, title = "Group A")
plot_graph(M_b, title = "Group B")
plot_graph(M_a, title = "", community = T)
plot_graph(M_b, title = "", community = T)

# Difference graph
A_means = colMeans(A_flat_boot)
B_means = colMeans(B_flat_boot)
G_delta = abs(A_means - B_means)

# Plot of differences at different thresholds
t_range = quantile(G_delta, probs = c(0.85, 0.95))
t_diff = seq(t_range[1], t_range[2], length.out = 4)
adj_matrix = function(t) as.integer(G_delta >= t)
delta_matrices = t(sapply(t_diff, adj_matrix))
titles = unlist(lapply(t_diff, function(t) paste(c("t =", round(t,2)), collapse = " ")))

plot_degree_distributions = function(flat_matrices, titles) {
  colors = c("deeppink", "deepskyblue", "gold2", "orchid")
  for(i in 1:4) {
    M = to_matrix(flat_matrices[i,])
    G = graph_from_adjacency_matrix(M, mode="undirected")
    deg = degree(G, mode="all")
    deg.dist = degree.distribution(G, cumulative=T, mode="all")
    if(i == 1) {
      plot(x=0:max(deg), y=1-deg.dist, pch=19, type="l", col=colors[i],
           xlab="Degree", ylab="Cumulative frequency",
           main="Degree distribution at different thresholds" , lwd=2)
    }
    else {
      lines(x=0:max(deg), y=1-deg.dist,
            col=colors[i], lwd=2)
    }
  }
  grid()
  legend("bottomright", legend=titles, col=colors, lty=1, lwd=2)
}

par(mfrow=c(1,1), mai=c(.5,.5,.5,.5))
plot_degree_distributions(delta_matrices, titles)