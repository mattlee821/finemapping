\name{run_susie_rss}
\alias{run_susie_rss}
\title{Run SuSiE-RSS}
\usage{
run_susie_rss(variants, R, effective_n = NULL, L = 10L,
  coverage = 0.95, min_abs_corr = 0.5, max_iter = 50L,
  refine = FALSE, psd_max_variants = 5000L)
}
\arguments{
\item{variants}{Validated variant data with signed \code{EZ} or valid
\code{ES}/\code{SE}, and \code{SS} unless \code{effective_n} is supplied.}
\item{R}{Signed LD matrix in the same allele-aware variant order.}
\item{effective_n}{Optional explicit scalar sample-size override. Otherwise
the minimum available per-variant \code{SS} is used.}
\item{L}{Maximum number of single-effect components.}
\item{coverage}{Credible-set posterior coverage.}
\item{min_abs_corr}{Minimum absolute correlation for reported credible sets,
or \code{NULL} to disable purity filtering.}
\item{max_iter}{Maximum number of IBSS iterations.}
\item{refine}{Whether to run the SuSiE refinement step.}
\item{psd_max_variants}{Maximum variant count for dense \eqn{O(M^3)} PSD
eigenvalue checking.}
}
\value{
A list containing variant PIPs, credible-set summaries and membership,
diagnostics, provenance, and the \code{coloc}-ready fitted \code{susieR}
model.
}
\description{
Validates summary statistics and LD, derives signed Z-scores only from
effect/standard-error inputs when needed, and fits the pinned \code{susieR}
RSS implementation. This method reports Bayesian PIPs and per-effect
credible sets.
}
