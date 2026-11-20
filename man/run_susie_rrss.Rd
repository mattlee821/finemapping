\name{run_susie_rrss}
\alias{run_susie_rrss}
\title{Run robust SuSiE-RSS}
\usage{
run_susie_rrss(variants, R, reference_sample_size, effective_n = NULL,
  L = 10L, coverage = 0.95, min_abs_corr = 0.5,
  max_iter = 100L, refine = FALSE, R_mismatch = "eb_mix",
  psd_max_variants = 5000L, eigen_max_variants = 20000L,
  allow_missing_n = FALSE)
}
\arguments{
\item{variants}{Validated variant data with signed \code{EZ} or valid
\code{ES}/\code{SE}, and \code{SS} unless \code{effective_n} is supplied or
\code{allow_missing_n = TRUE}.}
\item{R}{Signed LD matrix in the same allele-aware variant order.}
\item{reference_sample_size}{Positive LD reference sample size \eqn{B} used
for \code{R_finite}.}
\item{effective_n}{Optional explicit scalar sample-size override. Otherwise
the minimum available per-variant \code{SS} is used.}
\item{L}{Maximum number of single-effect components.}
\item{coverage}{Credible-set posterior coverage.}
\item{min_abs_corr}{Minimum absolute correlation for reported credible sets,
or \code{NULL} to disable purity filtering.}
\item{max_iter}{Maximum number of IBSS iterations.}
\item{refine}{Whether to run the SuSiE refinement step.}
\item{R_mismatch}{Robust LD-mismatch model: \code{"none"}, \code{"eb"}, or
\code{"eb_mix"}. The default is \code{"eb_mix"}.}
\item{psd_max_variants}{Maximum variant count for dense \eqn{O(M^3)} PSD
eigenvalue checking.}
\item{eigen_max_variants}{Maximum variant count for dense \eqn{O(M^3)}
\code{eigen(R)} evaluation of the \code{Q_art} subspace diagnostic; loci above
this size retain full \code{R_mismatch} (\code{lambda_bias} and SER
initialization) while recording \code{q_art_eigen = "skipped_large_m"}.}
\item{allow_missing_n}{Allow \code{z}-only fitting (\code{n = NULL}) when
neither \code{SS} nor \code{effective_n} is available.}
}
\value{
A list containing variant PIPs, credible-set summaries and membership,
robust-model diagnostics, provenance, and the \code{coloc}-ready fitted
\code{susieR} model.
}
\description{
Fits the robust RSS extension in \code{susieR} with
\code{R_finite = reference_sample_size} and configurable
\code{R_mismatch} (default \code{"eb_mix"}).
The \code{SER_fallback} reliability flag is recorded as a diagnostic; the
wrapper does not silently replace the full fit with the one-effect model.
}
