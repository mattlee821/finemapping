\name{read_ld_matrix}
\alias{read_ld_matrix}
\title{Read a square LD correlation matrix}
\usage{
read_ld_matrix(matrix_path, ids_path = NULL, scale = 1)
}
\arguments{
\item{matrix_path}{Path to a binary float64 NumPy \code{.npy} matrix or a
tab-delimited \code{.tsv} / \code{.tsv.gz} matrix.}
\item{ids_path}{Optional path to a one-identifier-per-line text file paired
with a \code{.npy} matrix. Defaults to \code{<stem>.ids.txt}.}
\item{scale}{Optional positive multiplier applied in-place while loading
(e.g. \code{n - 1} to construct \eqn{X^\top X = (n - 1)R} without allocating a
second \eqn{M \times M} matrix). Defaults to \code{1}.}
}
\value{
A numeric square matrix with row and column names set to the variant
identifiers.
}
\description{
Loads an aligned regional LD correlation matrix from either a binary
little-endian float64 \code{.npy} file plus companion identifier list or a
tab-delimited text file.
}
