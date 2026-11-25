library(MixedDeepMM)

# A small synthetic demonstration of the shared-tail layer.
sim <- simulate_mixeddeepmm(n = 40L, seed = 42L)
fit <- mixeddeepmm(sim$data, sim$types,
  family = "gaussian", heads = "two",
  layers = list(k = c(2L, 2L, 2L),
                r = list(continuous = 2L, discrete = 2L, tail = 1L)),
  control = list(seed = 42L, nstart = 1L, max_iter = 3L,
                 min_iter = 2L, mc_samples = 8L, mc_eval = 12L))
post <- predict(fit, type = "layer", cluster_level = "shared_tail")
clusters <- max.col(post, ties.method = "first")
print(table(clusters))
