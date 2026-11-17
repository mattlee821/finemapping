# Proposed architecture and data contracts

**Status:** The R package implements P-value ranking, p-LD, SuSiE-RSS, SuSiE-rRSS, GWAS-statistic validation, and LD-matrix validation. Synthetic tests and real-data smoke tests cover GWAS-VCF parsing and P-value ranking on three studies. Public-S3 Pan-UKB extraction produced a validated 700-variant EUR matrix for a 100-kb region. A bounded singleton Nextflow fine-mapping run completed on Ferkingstad P43320/CRYBB2 for three selected genome-wide-significant loci with LD-matched variants, and the gene-region output labelled 20,280 variants as cis. Method outputs are aggregated per GWAS as method-specific `.rda` files, while region-level TSV/JSON audit sidecars remain. The runner now defaults to all mapped clusters and queries each full positional cluster plus a flank; the historical three-locus test used lead-centered windows. At the strongest mapped locus, SuSiE-RSS still did not converge after 200 iterations. SuSiE-rRSS converged at iteration 124 with a 200-iteration cap, but its LD/summary-statistics reliability warning remained; `SER_fallback` is a flag, not replacement of the multi-effect fit. These diagnostics mean the results are not scientifically validated. The significance, clustering-distance, interval-flank, and cis-window parameters are provisional; full-genome and directory-scale throughput remain untested. The fine-mapping workflow currently accepts one GWAS at a time. See [the Pan-UKB LD guide](panukb-ld.md).

## System boundaries

GWASprocessing is the upstream producer of standardized GWAS-VCF files. GWASfinemapping consumes those files; it must not make GWASprocessing a runtime dependency or repeat its harmonization. This repository owns its validation, reference matching, methods, workflow, outputs, dependencies, and versioning.

The same R functions should be usable interactively and from scripts. Nextflow orchestrates tasks; it does not contain statistical method logic. Python is reserved for supporting components such as a validated `ldcov` or later PolyFun integration.

## Proposed components

| Component | Responsibility |
|---|---|
| R input/QC | Read VCF headers and records; validate required INFO/FORMAT fields, statistics, uniqueness, genome build, and metadata; preserve input fields and exclusion reasons. |
| R loci/regions | Define genome-wide loci or consume explicit intervals; gene annotation and cis/trans classification remain separate from fine-mapping. |
| R LD interface | Request regional Pan-UKB LD; match variants by coordinates and alleles; align signed correlations and Z-scores; validate and cache a uniquely keyed matrix. |
| R methods | Implement ordered Phase 1 methods and normalize shared result fields while retaining method-specific diagnostics. |
| Nextflow | Orchestrate GWAS × region × method work, reuse LD across compatible methods, allocate SLURM resources, and aggregate outputs/provenance. |
| Python support | Wrap established tools only where needed; no duplicate implementation of statistical methods. |

## Versioned run configuration proposal

Use one YAML run configuration with a required `schema_version`. Validate it against a checked-in JSON Schema before launching work. The proposed top-level contract is:

```yaml
schema_version: 1
run:
  id: example-run
  genome_build: GRCh38
gwas:
  id: trait-or-study-id
  vcf: /path/to/standardized.gwas.vcf.gz
  trait: example-trait
  ancestry: recorded-label
regions:
  mode: genome_wide # or explicit / gene_regions
  significance_threshold: null
  merge_distance_bp: null
  extension_bp: null
ld:
  reference: panukb
  release: required-release-id
  population: required-explicit-label
  reference_sample_size: null # required when using SuSiE-rRSS
  variant_index: /path/to/matched.variant_index.parquet
methods:
  - id: p_value
    top_n: null
  - id: p_ld
    r2_thresholds: [0.1, 0.01, 0.001]
sample_size:
  estimator: minimum_variant_ss
  effective_n_override: null
resources:
  cpus: 1
  memory: 4.GB
  time: 1.h
```

This is the version 1 interface proposal, mirrored in `config/schema.json` and `config/example.yaml`. Defaults currently confirmed by the project owner:

- p-LD\* thresholds: `r2 = 0.1, 0.01, 0.001`, each analyzed independently.
- Gene-region default: gene body extended by 1 Mb upstream and downstream.
- SuSiE sample-size estimator: minimum QC-passing per-variant `SS` in the region, unless an explicit effective-N override is supplied.
- `susieR`: pin GitHub commit `4134b03dceb0e7c0e1a8885a926781ba4b0d6743` (tag 2.0.0) for both RSS methods. SuSiE-rRSS uses `R_finite = B`; `R_mismatch` defaults to `"eb_mix"` and is configurable among the pinned implementation's supported choices.

Required validation rules:

- `schema_version` must be supported; unknown keys and invalid method settings must fail preflight with clear errors.
- The JSON Schema and example config are checked in; workflow-level preflight enforcement is not implemented yet.
- GWAS genome build and selected LD release/population are explicit. Never infer a Pan-UKB population from a GWAS ancestry label.
- Trait-specific overrides are explicit and recorded; they do not silently replace run-wide defaults.
- Paths, input checksums, reference release, population, sample size, parameters, software versions, and workflow revision are recorded in run provenance.
- A method-specific input or configuration failure must not be reported as a successful empty result.

## GWAS-VCF input contract proposal

The expected input is the standardized GWAS-VCF produced upstream, on GRCh38, with ALT as the effect allele. Preserve the original record and relevant metadata, including `ID`, `INFO/SNPID`, `INFO/HM`, `FORMAT/AFS`, and other present upstream annotations.

At minimum validate and normalize:

| Field | Contract |
|---|---|
| `CHROM`, `POS`, `REF`, `ALT` | Required variant identity and GRCh38 location; reject invalid coordinates and unsupported multiallelic records unless explicitly handled. |
| `FORMAT/ES`, `SE` | Effect estimate relative to ALT and its standard error when present; validate finite values and positive `SE`. |
| `FORMAT/LP` | `-log10(P)` when present; validate range and derive P-values stably only when needed. |
| `FORMAT/EZ` | Z-score when present; otherwise derive only from valid effect/SE data. |
| `FORMAT/SS` | Preserve per-variant sample size; do not silently collapse heterogeneous values to a scalar. |
| `FORMAT/AF`, `AFS`, `SI`, `NC` | Preserve and validate when present; do not require fields irrelevant to a particular study design. |

Preflight should detect duplicate variants, missing/invalid required statistics, inconsistent field declarations, and malformed records. It must state whether a statistic was supplied or derived. The precise required-field matrix by method, multi-sample VCF behavior, and acceptable sample-size heterogeneity remain to be tested against representative GWASprocessing outputs before finalizing the schema.

## Variant and LD contract

GWAS-VCF variants are GRCh38 while the specified Pan-UKB index is GRCh37. Regional LD requests must use allele-aware build conversion and retain a mapping record for every input variant, including mapped coordinates/alleles, transformations, status, and exclusion reason. Matching by rsID alone is prohibited.

Before any method consumes LD, validate dimensions, order, variant coverage, allele orientation, symmetry, diagonal, correlation range, missingness, and numerical stability/PSD. Any regularization or matrix repair must be explicit, configurable, and recorded. A matrix value outside the reference's stored band is unavailable LD, not zero LD. Unavailable, unmatched, and measured near-zero values must remain distinct.

Use a cache key containing at least the reference release, population, build mapping/version, region, extraction settings, variant set, and variant order. Validate cached contents before reuse. Align Z and R together; when genotype coding is reversed, use the corresponding row/column sign flips in R and sign change in Z.

## Proposed normalized result contract

Every per-region/per-method result should include common identifiers:

`schema_version`, `run_id`, `gwas_id`, `trait`, `region_id`, `method_id`, `method_version`, `CHROM`, `POS`, `REF`, `ALT`, `ID`, `SNPID`, `genome_build`, original GWAS statistics, `ld_reference`, `ld_population`, execution status, and provenance reference.

Method-specific fields are additive:

- Prioritization: rank, selected flag, selection threshold, and pruning/exclusion reason. P-value methods must not emit `PIP` or credible-set membership.
- Bayesian methods: PIP, credible-set/component identifier, set coverage, purity, convergence, and method-specific diagnostics.

Produce per-region/per-method results, aggregate genome-wide results, method diagnostics, explicit failure/exclusion reports, and run-level provenance. Use typed missing values for fields not defined by a method; do not encode “not applicable” as zero.

## Confirmed decisions and remaining evidence gaps

Confirmed and implemented:

1. p-LD uses independent `r^2 = 0.1, 0.01, 0.001` runs; top-N applies only to P-value prioritization.
2. `susieR` is pinned to GitHub tag 2.0.0's commit; synthetic tests exercise both RSS methods and diagnostics.
3. RSS uses the minimum non-missing regional `SS` unless an explicit effective-N override is supplied; per-variant values and the selected estimator are retained.

Remaining evidence gaps and design decisions:

1. Validate genome-wide locus thresholds/merge/extension, gene annotation source/build, and gene-region/cis-window defaults on target use cases.
2. Select and validate the allele-aware GRCh38→GRCh37 mapping source and ambiguous/indel/palindromic variant policy.
3. Validate binary-trait interpretation of `SS`/`NC` before treating minimum `SS` as SuSiE `n`; it remains an explicitly recorded approximation.
7. Confirm whether the initial release should error on incomplete regional LD or permit a documented alternative such as splitting/redefining the region.

