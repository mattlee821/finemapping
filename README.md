# GWASfinemapping

GWASfinemapping is an independent downstream project for statistical fine-mapping and variant prioritization from standardized GWAS-VCF summary statistics. It does not depend on or repeat GWASprocessing; GWAS-VCF files are consumed through the input contract described in [the architecture proposal](docs/architecture.md).

## Architecture

- Reusable statistical and validation functions in an R package.
- Nextflow workflow orchestration for one-GWAS runs and a separate directory preflight pilot, with local/HPC execution including SLURM.
- Python only for supporting components whose established implementations require it.
- Pan-UKB regional LD extraction with explicit reference population, release, variant matching and allele alignment.
- Shared, validated LD inputs across compatible methods, with provenance and diagnostics.

## Tool setup

Run `bash scripts/source.sh` to install the pinned micromamba environment and restore
the R dependencies. Then use `source scripts/source.sh` to activate it in a shell.
The environment and its pinned Python tools are built and published to GHCR as
`ghcr.io/mattlee821/finemapping`; see [the Pan-UKB LD guide](docs/panukb-ld.md)
for public-S3 index building and regional extraction.

## Development status

The literature review and design constraints are documented in [docs/literature_review.md](docs/literature_review.md) and [docs/architecture.md](docs/architecture.md). The R package implements P-value prioritization, p-LD, SuSiE-RSS, SuSiE-rRSS, and input/LD validation. Synthetic tests and real-data smoke tests cover GWAS-VCF parsing and P-value ranking on Locke, O'Mara, and Ferkingstad inputs. Public-S3 Pan-UKB EUR LD extraction was validated, followed by bounded end-to-end Nextflow runs on Ferkingstad P43320/CRYBB2: one locus and three selected genome-wide-significant loci with LD-matched variants. Method results are aggregated per GWAS in `methods/*.rda`, with TSV/JSON audit sidecars retained by region. The runner now defaults to all mapped significant-variant clusters and uses each cluster's full span plus a flank; the earlier real-data test was explicitly capped at three and used lead-centered windows. The CRYBB2 gene-region output labels 20,280 variants as cis. At the strongest mapped locus, SuSiE-RSS still did not converge after 200 iterations; SuSiE-rRSS converged at iteration 124 with a 200-iteration cap but retained its LD/summary-statistics reliability warning. Its `SER_fallback` field is a flag, not a replacement of the multi-effect fit. These runs validate workflow integration, not scientific reliability. The test settings for significance, cluster distance, interval flank, and cis window are not established production choices, and directory-scale throughput remains untested. See [docs/panukb-ld.md](docs/panukb-ld.md) for output structure, interval rules, and run limitations. The previous scripts, data, and README are preserved under [archive/legacy/](archive/legacy/).

Method order:

1. P-value prioritization.
2. P-value plus LD (`r^2`) prioritization.
3. SuSiE-RSS.
4. SuSiE-rRSS.
5. CARMA and PolyFun + SuSiE only after Phase 1 is validated.

P-value-based methods are prioritization only; they do not produce PIPs or Bayesian credible sets.
