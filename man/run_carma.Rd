\name{run_carma}
\alias{run_carma}
\title{Run CARMA Bayesian fine-mapping with outlier detection}
\usage{
run_carma(variants, R, lambda = 1,
  effect_size_prior = c("Spike-slab", "Cauchy"), tau = 0.04, L = 10L,
  coverage = 0.95, min_abs_corr = 0.5, outlier_switch = TRUE,
  outlier_bf_index = 1 / 3.2, max_model_dim = 200000L, all_iter = 3L,
  all_inner_iter = 10L, epsilon_threshold = 1e-5,
  psd_max_variants = 5000L, seed = 1L)
}
\arguments{
\item{variants}{Validated variant data with signed \code{EZ} or valid
\code{ES}/\code{SE}.}
\item{R}{Signed LD correlation matrix in the same allele-aware variant order.}
\item{lambda}{Poisson model-size hyperparameter \eqn{\eta} (default \code{1}).}
\item{effect_size_prior}{Effect-size prior: \code{"Spike-slab"} (default) or
\code{"Cauchy"}.}
\item{tau}{Prior precision parameter for the Spike-slab effect-size prior
(default \code{0.04}).}
\item{L}{Maximum number of causal variants per locus (default \code{10}).}
\item{coverage}{Credible-set posterior coverage \eqn{\rho} (default
\code{0.95}).}
\item{min_abs_corr}{Minimum pairwise absolute correlation threshold for
credible sets (default \code{0.5}).}
\item{outlier_switch}{Whether to enable CARMA's Bayesian LD/Z-score discrepancy
outlier detection (default \code{TRUE}).}
\item{outlier_bf_index}{Bayes-factor threshold for flagging discrepancy
outliers (default \code{1 / 3.2}).}
\item{max_model_dim}{Maximum number of top candidate models retained.}
\item{all_iter}{Maximum outer Shotgun iterations.}
\item{all_inner_iter}{Maximum inner Shotgun epochs per outer iteration.}
\item{epsilon_threshold}{Per-variant convergence tolerance.}
\item{psd_max_variants}{Maximum variant count for dense \eqn{O(M^3)} PSD
checking.}
\item{seed}{Integer RNG seed for deterministic Shotgun stochastic search.}
}
\value{
A list containing variant PIPs and \code{CARMA_OUTLIER} flags, credible-set
summaries and membership, outlier/convergence diagnostics, and provenance.
}
\description{
Executes CARMA (Yang et al., 2023; \code{ZikunY/CARMA}) Shotgun stochastic
search with Poisson model-space priors and Bayesian hypothesis testing for
summary-statistic / external-LD discrepancies, avoiding dense sparse-matrix
coercions and redundant \eqn{M \times M} copies on large biobank loci.
}
