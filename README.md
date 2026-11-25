# MixedDeepMM 0.1.0

Deep Gaussian, Student-t and skew-t mixture models for continuous observations
and mixed-type measurements (continuous, count, binary, ordinal and nominal).
Continuous-only data use the integrated continuous model engines; mixed data
use a single-head or two-head Monte Carlo EM model. The two-head model connects
continuous and discrete observation heads through a shared mixture tail.

## Installation

Install the source archive from its local directory:

```r
install.packages("MixedDeepMM_0.1.0.tar.gz", repos = NULL, type = "source")
library(MixedDeepMM)
packageVersion("MixedDeepMM")
```

Optional packages are required only by the corresponding methods:
`mvtnorm` and `corpcor` for continuous DGMM, `GIGrvg` for skew-t fitting,
`cluster` for Gower dissimilarity/silhouette, and `mclust` for a supported
continuous-model initialisation choice. Install any that you need with
`install.packages()`.

## Continuous observations

```r
library(MixedDeepMM)
X <- as.data.frame(scale(iris[, 1:4]))
fit <- mixeddeepmm(X, family = "t", heads = "single",
  layers = list(k = c(3L, 2L), r = c(2L, 1L)),
  control = list(seed = 42L, nstart = 1L, nu = 8,
                 estimate_nu = TRUE, max_iter = 10L))
g <- predict(fit, type = "cluster", cluster_level = "first")
mixeddeepmm_ari(iris$Species, g)
```

## Mixed observations and shared-tail clustering

```r
sim <- simulate_mixeddeepmm(n = 40L, seed = 42L)
fit <- mixeddeepmm(sim$data, types = sim$types,
  family = "gaussian", heads = "two",
  layers = list(k = c(2L, 2L, 2L),
                r = list(continuous = 2L, discrete = 2L, tail = 1L)),
  control = list(seed = 42L, nstart = 1L, max_iter = 3L,
                 min_iter = 2L, mc_samples = 8L, mc_eval = 12L))
post <- predict(fit, type = "layer", cluster_level = "shared_tail")
clusters <- max.col(post, ties.method = "first")
D <- mixeddeepmm_gower_distance(sim$data, sim$types)
mixeddeepmm_gower_silhouette(clusters, D)
```

For a single-head mixed model, `cluster_level = "first"` selects the common
first mixture layer. For a two-head mixed model, `"shared_tail"` selects its
last shared mixture layer. The previous `"paper"` value remains an alias and
selects the appropriate layer based on the architecture.

## Statistical scope

The mixed-data estimator is an adapted, shallow conditional-measurement MCEM
model; it is not an exact reproduction of the published multi-layer MDGMM,
NSEP initialisation, pruning, or within-iteration silhouette selection.
Student-t and skew-t are separate distributional extensions. The Monte Carlo
likelihood and BIC are approximations. Gower silhouette and ARI/MR are
computed separately from fitting; true class labels must not be used when
selecting an unsupervised partition.

See `?mixeddeepmm`, `?predict.mixeddeepmm`, and
`inst/references/PROVENANCE.md` for definitions and source attribution.

## Licence and CRAN release

MixedDeepMM 0.1.0 is released under **GPL (>= 3)**. Code adapted from the
original GPL-licensed `deepgmm` project retains upstream author attribution.
See `inst/COPYRIGHTS` and `inst/references/PROVENANCE.md`.

For a CRAN submission, run `R CMD build MixedDeepMM` and
`R CMD check --as-cran MixedDeepMM_0.1.0.tar.gz` in a local R environment.
A CRAN submission also requires checking the package against the current
CRAN repository policy, including contributor rights and test results.
