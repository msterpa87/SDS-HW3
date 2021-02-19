library(DescTools)
library(igraph)
library(matrixStats)

load("hw3_data.RData")


average_cor = function(M) FisherZInv(colMeans(FisherZ(M)))
to_matrix = function(M_flat) matrix(as.integer(M_flat), nrow=116, byrow=T)

# test for empty intersection of [-t,t] with CI(j,k)
int.test <- function(lower, upper, t) (t < lower) | (-t > upper)

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
  A_matrix = get_adjancency_matrix(lower = A_fisher[2,], upper = A_fisher[3,],
                                   threshold = t_a)
  B_matrix = get_adjancency_matrix(lower = B_fisher[2,], upper = B_fisher[3,],
                                   threshold = t_b)
  
  return(list(A = A_matrix, B = B_matrix))
}

get_adjancency_matrix = function(lower, upper, threshold) {
  M = to_matrix(sapply(1:length(lower), function(i) int.test(lower[i], upper[i], threshold)))
  M[is.na(M)] = 0
  return(M)
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


b = 100

# point estimate correlation
cor_a = average_cor(A_flat)
cor_b = average_cor(B_flat)
cor_diff = cor_a - cor_b

t = mean(cor_diff, na.rm = T);t

A_boot = matrix(NA, nrow = b, ncol = 116^2)
B_boot = matrix(NA, nrow = b, ncol = 116^2)
L_boot_a = matrix(NA, nrow = b, ncol = 116^2)
U_boot_a = matrix(NA, nrow = b, ncol = 116^2)
L_boot_b = matrix(NA, nrow = b, ncol = 116^2)
U_boot_b = matrix(NA, nrow = b, ncol = 116^2)

n = 12

d = 116
m = d*(d-1)/2
alpha = 0.05/m
z = qnorm(1-alpha)

for(i in 1:b) {
  idx_a = sample(1:n, n, replace = T)
  idx_b = sample(1:n, n, replace = T)
  
  A_boot[i,] = average_cor(A_flat[idx_a,])
  B_boot[i,] = average_cor(B_flat[idx_b,])
  
  fisher_a = sapply(A_boot[i,], function(x) CorCI(x, n = 145, conf.level = 1-alpha))
  L_boot_a[i,] = fisher_a[2,]
  U_boot_a[i,] = fisher_a[3,]
  
  fisher_b = sapply(B_boot[i,], function(x) CorCI(x, n = 145, conf.level = 1-alpha))
  L_boot_b[i,] = fisher_b[2,]
  U_boot_b[i,] = fisher_b[3,]
}

D_boot = A_boot - B_boot

# Normal confidence intervals
se_diff = apply(D_boot, MARGIN = 2, sd)
se_diff[is.na(se_diff)] = 0

lower_a = colMeans(L_boot_a)
upper_a = colMeans(U_boot_a)
lower_b = colMeans(L_boot_b)
upper_b = colMeans(U_boot_b)

# Adjacency matrices
M_a = get_adjancency_matrix(lower_a, upper_a, t_a)
M_b = get_adjancency_matrix(lower_b, upper_b, t_b)

par(mfrow=c(2,2), mai=c(.1,.1,.2,.1))
plot_graph(M_a, title = "Group A")
plot_graph(M_b, title = "Group B")
plot_graph(M_a, title = "", community = T)
plot_graph(M_b, title = "", community = T)

D_boot = A_boot - B_boot
L_boot_diff = matrix(NA, nrow = b, ncol = 116^2)
U_boot_diff = matrix(NA, nrow = b, ncol = 116^2)

lower_diff = lower_a - lower_b
upper_diff = colMeans(U_boot_diff)

t_diff = mean(D_boot, na.rm = T)
M_diff = get_adjancency_matrix(lower_diff, upper_diff, t_diff)
par(mfrow=c(1,2), mai=c(.1,.1,.2,.1))
plot_graph(M_diff, title = "Difference Graph")
plot_graph(M_diff, title = "", community = T)
