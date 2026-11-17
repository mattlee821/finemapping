\name{validate_ld_matrix}
\alias{validate_ld_matrix}
\title{Validate an LD correlation matrix}
\usage{
validate_ld_matrix(R, variant_ids = NULL, check_psd = FALSE,
  tolerance = 1e-8, psd_max_variants = 5000L)
}
\arguments{
\item{R}{Numeric square correlation matrix with variant identifiers as row
and column names.}
\item{variant_ids}{Optional expected variant order.}
\item{check_psd}{Whether to check positive semidefiniteness.}
\item{tolerance}{Positive tolerance for numerical validation.}
\item{psd_max_variants}{Maximum matrix dimension for dense \eqn{O(M^3)}
eigenvalue verification when \code{check_psd = TRUE}. Matrices larger than
this threshold retain symmetry, diagonal, finite, and range validation while
recording \code{psd_status = "skipped_large_m"}.}
}
\value{
Invisibly returns \code{TRUE} (with attribute \code{psd_status}) if the
matrix passes validation; otherwise raises an error with the failed condition.
}
\description{
Checks dimensions, identifier order, missingness, symmetry, diagonal values,
correlation range, and optionally positive semidefiniteness. No matrix repair
is performed.
}
