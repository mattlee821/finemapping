\name{run_polyfun_susie}
\alias{run_polyfun_susie}
\title{Run functionally informed PolyFun + SuSiE fine-mapping}
\usage{
run_polyfun_susie(variants, R, prior_weights = NULL, polyfun_priors = NULL,
  missing_prior_policy = c("baseline_floor", "drop", "error"),
  effective_n = NULL, L = 10L, coverage = 0.95, min_abs_corr = 0.5,
  max_iter = 100L, refine = FALSE, psd_max_variants = 5000L,
  allow_missing_n = FALSE)
}
\arguments{
\item{variants}{Validated variant data with signed \code{EZ} or valid
\code{ES}/\code{SE}, GRCh37 \code{LD_CHROM}/\code{LD_POS}/\code{LD_REF}/\code{LD_ALT}
(or \code{CHROM}/\code{POS}/\code{REF}/\code{ALT}), and \code{SS} unless
\code{effective_n} is supplied or \code{allow_missing_n = TRUE}.}
\item{R}{Signed LD matrix in the same allele-aware variant order.}
\item{prior_weights}{Optional explicit non-negative numeric vector of prior
causal weights. When omitted, priors are resolved from \code{polyfun_priors} or
a \code{POLYFUN_SNPVAR} column.}
\item{polyfun_priors}{Optional path to an installed PolyFun reference directory
(\code{by_chrom/chr1..22.parquet}) or a custom PolyFun per-SNP heritability
table (\code{CHR, BP, A1, A2, SNPVAR}). Defaults to
\code{$GWASFM_REFERENCE_ROOT/polyfun}.}
\item{missing_prior_policy}{How to handle autosomal variants absent from the
PolyFun prior table: \code{"baseline_floor"} (default; assigns the PolyFun
\eqn{q = 100} floor so no candidate causal variant is dropped), \code{"drop"},
or \code{"error"}.}
\item{effective_n}{Optional explicit scalar sample-size override.}
\item{L}{Maximum number of single-effect components.}
\item{coverage}{Credible-set posterior coverage.}
\item{min_abs_corr}{Minimum absolute correlation for reported credible sets.}
\item{max_iter}{Maximum number of IBSS iterations.}
\item{refine}{Whether to run the SuSiE refinement step.}
\item{psd_max_variants}{Maximum variant count for dense \eqn{O(M^3)} PSD
eigenvalue checking.}
\item{allow_missing_n}{Allow \code{z}-only fitting (\code{n = NULL}) when
neither \code{SS} nor \code{effective_n} is available.}
}
\value{
A list containing functionally informed variant PIPs (plus \code{SNPVAR},
\code{PRIOR_WEIGHT}, and \code{PRIOR_MATCHED}), credible-set summaries and
membership, prior-coverage diagnostics, provenance, and the \code{coloc}-ready
fitted \code{susieR} model.
}
\description{
Implements PolyFun + SuSiE (Weissbrod et al., 2020) by aligning functional
per-SNP heritability priors (\code{SNPVAR} with \eqn{q = 100} max/min ratio
truncation) onto GRCh37 locus alleles and fitting \code{susieR} with
\code{prior_weights}.
}
