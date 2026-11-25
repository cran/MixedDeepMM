\name{mixeddeepmm}
\alias{mixeddeepmm}
\title{Fit deep mixture models for continuous or mixed responses}
\description{A unified model interface using existing continuous DGMM, RDMM and DStMM engines and a new single- or two-head mixed-response Monte Carlo EM estimator. Mixed Student-t and skew-t kernels extend the Gaussian formulation.}
\usage{
mixeddeepmm(
  data,
  types = NULL,
  family = c("gaussian", "t", "skew-t"),
  heads = c("single", "two"),
  layers,
  control = list()
)
}
\arguments{
\item{data}{Data frame or all-numeric matrix of observations. Missing values are not supported.}
\item{types}{Named character vector assigning every column one of \code{continuous}, \code{binary}, \code{ordinal}, \code{nominal}, or \code{count}. Defaults to all continuous for an all-numeric dataset.}
\item{family}{Latent kernel: Gaussian, Student-t or generalised-hyperbolic skew-t.}
\item{heads}{One shared latent head, or separate continuous/discrete heads with shared deeper tail.}
\item{layers}{List with \code{k} component counts and \code{r} dimensions. Single-head mixed example: \code{list(k=c(2,2),r=c(2,1))}; two-head example: \code{list(k=c(2,2,2),r=list(continuous=2,discrete=2,tail=1))}.}
\item{control}{Fitting controls: \code{max_iter} (35), \code{min_iter} (8), \code{tol} (0.002), \code{mc_samples} (48), \code{mc_eval} (192), \code{nstart} (1), \code{seed} (123), \code{nu} (8), \code{estimate_nu} (TRUE), \code{psi_floor} (0.03), \code{ridge} (0.001), \code{verbose} (FALSE), \code{init} ("kmeans"), \code{init_est} ("factanal"). Aliases \code{it}, \code{eps} and \code{M} are accepted. Mixed-response fits always use deterministic data-driven initialization; \code{seed} controls MCEM draws and the choice of deterministic initial partition, and \code{nstart} explores alternative deterministic variants. \code{init} and \code{init_est} apply only to legacy continuous-only engines.}
}
\details{For mixed responses the latent Gaussian hierarchy and measurement model are fitted using posterior draws conditional on continuous responses, importance weighted by non-Gaussian response likelihoods. Each fit saves a separate Monte Carlo approximation to the observed log likelihood for model comparison, with a delta-method integration-error estimate. The mixed Gaussian estimator uses the original paper's conditional Gaussian and GLLVM design with a component-specific two-head continuous measurement and one shallow Gaussian-mixture tail; it is not a full reimplementation of the published MDGMM estimation algorithm. Its structured PCA/indicator, clustering, and regression-based initialization is not the published NSEP algorithm; it uses no random initial parameter draws. MCEM integration still involves random sampling. Numerical integration and estimation sensitivity should be examined for the mixed Student-t and skew-t extensions. BIC counts the represented coefficients under reference-location constraints; residual factor rotational non-identifiability may affect interpretation.}
\value{A fitted object of class \code{mixeddeepmm}, supporting \code{print}, \code{summary}, \code{predict}, \code{logLik}, \code{fitted}, and \code{plot}.}
\examples{
sim <- simulate_mixeddeepmm(n = 40, seed = 1)
fit <- mixeddeepmm(sim$data, sim$types, layers = list(k = c(2,2), r = c(2,1)),
  control = list(max_iter = 2, min_iter = 2, mc_samples = 8, mc_eval = 16))
print(fit)
}
