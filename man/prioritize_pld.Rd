\name{prioritize_pld}
\alias{prioritize_pld}
\title{Prioritize variants using P-value-ordered LD pruning}
\usage{
prioritize_pld(variants, R, r2_thresholds)
}
\arguments{
\item{variants}{Data frame containing allele-aware variant fields and a P-value,
\code{LP}, or signed Z-score.}
\item{R}{Signed LD correlation matrix in exactly the same variant order.}
\item{r2_thresholds}{One or more unique values in \eqn{(0, 1]} specifying
\eqn{r^2} thresholds.}
}
\value{
A data frame with one row per variant per threshold, rank, selection flag,
threshold, and the selected variant responsible for each pruning decision.
}
\description{
For each threshold independently, visits variants in P-value order and selects
a variant only when its squared correlation with every previously selected
variant is below the threshold. Missing LD is an error, not zero LD.
}
\details{
This greedy procedure is a prioritization heuristic. It does not produce
posterior inclusion probabilities or Bayesian credible sets.
}
