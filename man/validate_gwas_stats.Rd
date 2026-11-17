\name{validate_gwas_stats}
\alias{validate_gwas_stats}
\title{Validate and normalize GWAS variant statistics}
\usage{
validate_gwas_stats(variants, require_p = FALSE, require_z = FALSE)
}
\arguments{
\item{variants}{Data frame with \code{CHROM}, \code{POS}, \code{REF}, and
\code{ALT}, and available GWAS statistics.}
\item{require_p}{Whether every row must have a finite log P-value.}
\item{require_z}{Whether every row must have a signed Z-score.}
}
\value{
A validated data frame. Missing \code{LP} values can be derived from P-values
or signed Z-scores; missing \code{EZ} values can be derived from \code{ES}/\code{SE}.
}
\description{
Checks allele-aware variant uniqueness, coordinates, alleles, and available
numeric GWAS statistics. It never infers a signed Z-score from an unsigned
P-value.
}
