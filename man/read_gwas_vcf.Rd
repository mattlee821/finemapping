\name{read_gwas_vcf}
\alias{read_gwas_vcf}
\title{Read a standardized GWAS-VCF}
\usage{
read_gwas_vcf(path)
}
\arguments{
\item{path}{Path to a VCF or gzipped VCF containing one GWAS sample.}
}
\value{
A \code{gwas_vcf} list containing the parsed \code{variants} data frame,
header metadata, sample identifier, and normalized input path.
}
\description{
Reads standard VCF fixed fields, INFO metadata, and sample FORMAT fields.
Numeric GWAS fields are validated, and the original INFO and FORMAT values
are preserved.
}
\details{
Exactly one sample column is supported per file. Multiallelic variants are
rejected by \code{\link{validate_gwas_stats}} rather than split or silently
discarded. A FORMAT tag that collides with a fixed VCF field or INFO-derived
\code{SNPID}/\code{HM} column is retained with a \code{FORMAT_} prefix; for
example, FORMAT/ID is returned as \code{FORMAT_ID}.
}
