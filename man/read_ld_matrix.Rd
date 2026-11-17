\name{read_ld_matrix}
\alias{read_ld_matrix}
\title{Read a square LD correlation matrix}
\usage{
read_ld_matrix(matrix_path, ids_path = NULL)
}
\arguments{
\item{matrix_path}{Path to a binary float64 NumPy \code{.npy} matrix or a
tab-delimited \code{.tsv} / \code{.tsv.gz} matrix.}
\item{ids_path}{Optional path to a one-identifier-per-line text file paired
with a \code{.npy} matrix. Defaults to \code{<stem>.ids.txt}.}
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
