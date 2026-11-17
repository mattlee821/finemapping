\name{prioritize_p}
\alias{prioritize_p}
\title{Rank variants by GWAS P-value}
\usage{
prioritize_p(variants, top_n = NULL)
}
\arguments{
\item{variants}{Data frame containing allele-aware variant fields and a P-value,
\code{LP}, or signed Z-score.}
\item{top_n}{Optional positive integer limiting the selected variants.}
}
\value{
A data frame with original and normalized statistics, deterministic rank,
selection flag, and selection reason.
}
\description{
Ranks variants within the supplied region by increasing P-value, using stable
log-P calculations when available. This is prioritization, not Bayesian
fine-mapping; no PIPs or credible sets are produced.
}
