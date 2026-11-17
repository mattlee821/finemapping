# Fine-mapping literature review

**Status:** Literature findings are complete for the locally supplied full texts and accessible primary sources. Several linked supplementary files and the CARMA article full text could not be retrieved; those gaps are recorded below. Recommendations are provisional where they depend on those gaps or on project choices.

## Scope and evidence

The review covers the five papers supplied in `docs/papers/`, the primary SuSiE-RSS and PolyFun papers, current `susieR` documentation, CARMA implementation documentation, RSparsePro as a relevant LD-mismatch comparator, and GeneHackman and `ldcov` as implementation references. Literature describes statistical behavior; the Pan-UKB release facts and GWAS-VCF field contract below are project requirements and still need release-specific preflight verification.

### Supplied-paper inventory

| File | Paper |
|---|---|
| [bbag574.pdf](./papers/bbag574.pdf) | Chen et al. (2026), “Systematic benchmarking of multi-ancestry fine-mapping under variant heterogeneity and LD mismatch,” *Briefings in Bioinformatics*, DOI [10.1093/bib/bbag574](https://doi.org/10.1093/bib/bbag574). |
| [ddaa148.pdf](./papers/ddaa148.pdf) | Hutchinson et al. (2020), “Fine-mapping genetic associations,” *Human Molecular Genetics*, DOI [10.1093/hmg/ddaa148](https://doi.org/10.1093/hmg/ddaa148). |
| [main.pdf](./papers/main.pdf) | Benner et al. (2017), “Prospects of Fine-Mapping Trait-Associated Genomic Regions by Using Summary Statistics from GWAS,” *The American Journal of Human Genetics* 101:539–551, DOI [10.1016/j.ajhg.2017.08.012](https://doi.org/10.1016/j.ajhg.2017.08.012). |
| [nihms-2112025.pdf](./papers/nihms-2112025.pdf) | Li and Zhou (2025), “Towards improved fine-mapping of candidate causal variants,” *Nature Reviews Genetics*, DOI [10.1038/s41576-025-00869-4](https://doi.org/10.1038/s41576-025-00869-4). |
| [nihms976173.pdf](./papers/nihms976173.pdf) | Schaid et al. (2018), “From genome-wide associations to candidate causal variants by statistical fine-mapping,” *Nature Reviews Genetics* 19:491–504, DOI [10.1038/s41576-018-0016-z](https://doi.org/10.1038/s41576-018-0016-z). |

The complete text of each supplied PDF was reviewed. The supplementary material is not all present in the supplied directory. PMC lists two Benner supplement PDFs and a Schaid DOCX, but their download endpoints returned access interstitials during this review; their contents are therefore **not** represented as reviewed. Chen et al.'s supplementary availability was not verified. These are literature-review limitations, not evidence that supplements do not exist.

## Findings relevant to GWASfinemapping

### Fine-mapping is not association ranking

A small P-value identifies association, not a causal variant. LD can make many variants similarly associated, and the lead association need not be the causal variant. Bayesian fine-mapping combines association evidence, an LD model, and prior assumptions to estimate uncertainty over candidate causal variants. PIPs and credible sets are conditional on that model and its inputs; they do not prove biological causality.

The supplied reviews emphasize that the meaning of a credible set depends on its construction. SuSiE's credible sets are per-effect sets intended to contain the causal variant for a particular single-effect component at a chosen posterior coverage, subject to the model and filtering rules. They should not be described as a single set covering all causal variants in a locus. Report the coverage and purity rules, and preserve the component association and diagnostics. P-value-only and p-value/LD-pruning results are prioritization outputs and must not be presented as PIPs or Bayesian credible sets.

### LD reference quality and mismatch

Summary-statistic methods need LD consistent with the GWAS statistics: variant identity, order, alleles, ancestry, and sample context all matter. Benner et al. show that reference-panel requirements grow with GWAS sample size. Their simulations found that roughly 1,000 target-population reference samples could be adequate for GWAS sample sizes up to about 10,000, while substantially larger panels were needed for larger GWAS; small generic panels such as 1000 Genomes can materially degrade fine-mapping. This is not a universal sample-size threshold for modern methods or Pan-UKB; it is evidence that reference size and compatibility must be recorded and tested.

The supplied reviews and SuSiE-RSS paper describe consequences of LD/statistic inconsistency, including distorted posterior probabilities, unreliable credible-set coverage, and spurious signals. Standard SuSiE-RSS diagnostics can detect some inconsistencies but do not make an arbitrary reference panel equivalent to in-sample LD. Allele flips are particularly dangerous: a sign mismatch in Z must be paired with the corresponding row and column sign changes in R. Matching by rsID alone is insufficient.

The 2026 Chen et al. benchmark reports that fine-mapping performance under variant heterogeneity and LD mismatch varies by method and scenario. Joint multi-ancestry analyses can gain information for shared signals but may lose ancestry-specific variants when variants are not shared across input sets; post-hoc per-ancestry analyses can retain those signals. The benchmark supports treating ancestry and variant intersection as explicit analysis choices, not assuming one ancestry-specific LD panel is appropriate for every meta-analysis.

Zhang et al.'s requested paper is identifiable as “Robust fine-mapping in the presence of linkage disequilibrium mismatch,” *Nature Genetics* 58 (2026), DOI [10.1038/s41588-026-02780-y](https://doi.org/10.1038/s41588-026-02780-y); its 2024 bioRxiv preprint is [10.1101/2024.10.29.620968](https://doi.org/10.1101/2024.10.29.620968). The version-of-record full text and preprint supplement ZIP were not accessible in this review; findings were checked against the preprint and the [official RSparsePro documentation](https://github.com/zhwm/RSparsePro_LD). RSparsePro models latent Z-scores to reduce LD inconsistency and reports PIPs and credible-set membership. It is useful evidence that mismatch remains a distinct modeling problem, but is **not** in the requested implementation scope.

### SuSiE-RSS and robust SuSiE-RSS

Zou et al. (2022), [“Fine-mapping from summary data with the ‘Sum of Single Effects’ model”](https://doi.org/10.1371/journal.pgen.1010299), describes the summary-statistic extension. SuSiE represents the vector of variant effects as a sum of up to `L` single-effect vectors; iterative Bayesian stepwise selection (IBSS) fits the model. This permits multiple causal effects in a region, unlike a single-causal-variant model. RSS uses summary statistics such as Z-scores and an LD matrix as an approximation to the individual-level sufficient-statistic likelihood. With exact in-sample sufficient statistics it can reproduce the corresponding individual-level analysis; external LD makes this an approximation. The model was developed primarily for quantitative-trait regression; case-control and meta-analysis inputs require careful interpretation.

The primary output includes per-variant PIPs, per-effect credible sets, and fit/convergence diagnostics. `L` is an upper bound, not a known causal count: too small a value can constrain the model; too large a value can produce weak or redundant components. The project should retain convergence information, credible-set purity, component-level details, and sensitivity diagnostics rather than treating every returned set as equally informative. RSS computation is designed to work from summary data without requiring individual-level genotypes; run time and memory still depend on region size, `L`, and the numerical work required by the fit.

Current `susieR` GitHub documentation (tag 2.0.0) describes `R_finite` and `R_mismatch`; the checked CRAN 0.14.2 documentation and GitHub 0.16.6 documentation did not contain these options. In this interface:

- `R_finite` can account for estimation uncertainty from a finite LD reference. It accepts a numeric reference-panel size; `TRUE` attempts to infer the size from `X`; and `FALSE` disables the adjustment, intended for in-sample LD.
- `R_mismatch = "eb_mix"` estimates regional mismatch and uses a two-component residual mixture to estimate variant-level mismatch probabilities.
- `SER_fallback` is an output/diagnostic alias for `R_reliability_flag`, not a switch that automatically reruns or replaces the model with SER. A separate one-effect SER rerun, if desired, must be explicit and separately recorded.

Thus, robust SuSiE-RSS is implemented as a distinct method using the same `susie_rss` entry point, pinned to commit `4134b03dceb0e7c0e1a8885a926781ba4b0d6743` (tag 2.0.0). Deterministic synthetic tests confirm the robust arguments and diagnostics are available in that version; representative-locus validation remains outstanding. The separate publication identity of the “SuSiE-rRSS paper” was not verified; the current mismatch vignette names SuSiE-rRSS as the robust extension. The project should describe it as a distinct configured method and provenance record without asserting an unverified separate paper.

### CARMA

Yang et al. (2023), “CARMA is a new Bayesian model for fine-mapping in genome-wide association meta-analyses,” *Nature Genetics* 55:1057–1065, DOI [10.1038/s41588-023-01392-0](https://doi.org/10.1038/s41588-023-01392-0), is the requested primary paper. The full article was not accessible from the checked publisher/Europe PMC sources, so a complete account of its model and computational properties is not established here. The [official CARMA package manual](https://github.com/ZikunY/CARMA/blob/master/man/CARMA.Rd) and [README](https://github.com/ZikunY/CARMA) document Z-statistics, signed LD correlations, optional functional annotations, variant-level PIPs, credible-set/model output, and LD-outlier identification; they also document user parameters including a default maximum of 10 causal variants and an outlier Bayes-factor threshold. These implementation claims must not be mistaken for a full review of the paper's derivations or limitations.

CARMA remains a Phase 2 candidate, not an implementation decision. Before integration, obtain and review the full paper and any available supplement, then test the established implementation on representative data. Do not infer that outlier handling resolves all forms of population, allele, variant, or reference-size mismatch.

### PolyFun + SuSiE

Weissbrod et al. (2020), “Functionally informed fine-mapping and polygenic localization of complex trait heritability,” *Nature Genetics* 52:1355–1363, DOI [10.1038/s41588-020-00735-5](https://doi.org/10.1038/s41588-020-00735-5), introduces PolyFun as a genome-wide prior-estimation framework, not a per-locus fine-mapper. It uses genome-wide functional annotations and regularized stratified LD-score regression, cross-fitting and heritability-based SNP bins to estimate per-SNP heritability and derive variant-level priors. A method such as SuSiE then consumes the priors for locus-level fine-mapping.

PolyFun therefore requires genome-wide summary/annotation/reference compatibility, not just a list of variants at a locus. Annotation genome build, variant coverage, LD-score reference, trait definition, prior estimation, and compatibility with the SuSiE variant set must be tracked. The original paper's principal analyses used in-sample LD; external-LD mismatch remains a concern. The [official PolyFun repository](https://github.com/omerwe/polyfun) currently describes a Python workflow with R `susier` and optional other fine-mappers. Its README recommends mutable `master` and notes that its lockfile may be outdated, so dependency pinning and software validation are prerequisites. Its separate supplements were linked but not extracted here. This is also Phase 2 only.

## Method comparison

The first four methods are the requested Phase 1 sequence. CARMA and PolyFun + SuSiE are Phase 2 and must wait until Phase 1 is validated.

| Method | Model / inputs | LD assumptions and mismatch | Outputs | Compute profile and principal limitations |
|---|---|---|---|---|
| **P-value prioritization** | Rank valid GWAS-VCF variants within a region using P or `LP`; optional top-N selection. | No LD required. | Original GWAS statistics, deterministic rank and selection status. **No PIPs or credible sets.** | Lightweight sorting. Does not distinguish correlated candidates or quantify causal uncertainty. |
| **p-LD\*** | Greedy P-value-ordered pruning using pairwise \(r^2\); input statistics plus regional LD. | Requires usable, correctly matched \(r^2\) for the comparisons made. Missing/out-of-band LD is unknown, not zero. Sensitive to threshold and variant set. | Selected variants per \(r^2\) threshold, rank and pruning/exclusion reason. **No PIPs or credible sets.** | Lightweight to moderate, depending on region and LD representation. Heuristic; thresholds are not posterior uncertainty and results can vary with tie-breaking and region boundaries. |
| **SuSiE-RSS** | Bayesian sum-of-single-effects regression; aligned Z, LD matrix, sample size, `L` and priors/settings. | External LD approximates in-sample LD; requires compatible variants, allele coding, ancestry and sufficient accuracy. Mismatch can distort PIPs and coverage. | PIPs, per-effect credible sets, purity and model/convergence diagnostics. | Summary-data scalable; cost grows with region size and `L`. Model and credible sets depend on correct inputs and specified priors. |
| **SuSiE-rRSS** | SuSiE-RSS with finite-reference and LD-mismatch modeling in pinned `susieR` 2.0.0. | Uses the same aligned Z/R requirements; `R_finite` and `R_mismatch` address particular finite-panel and mismatch uncertainty, not arbitrary absent/out-of-band LD. | SuSiE outputs plus robust settings, mismatch/reliability diagnostics and `SER_fallback`/`R_reliability_flag`. | Similar family of computation to RSS with additional estimation. Synthetic API tests pass; representative-locus comparison is outstanding. The flag is diagnostic, not automatic fallback. |
| **CARMA** (Phase 2) | Bayesian summary-statistic method; Z, signed LD; optional annotations per implementation docs. | Explicit LD-outlier functionality documented; full paper review is still needed to characterize assumptions and scope. | PIPs, credible/model output and outlier results per manual. | Not yet characterized sufficiently from primary text to estimate reliably; do not integrate based only on package manual. |
| **PolyFun + SuSiE** (Phase 2) | Genome-wide annotation-based prior estimation followed by locus-level SuSiE. | Requires consistent genome-wide annotations, LD scores/reference and fine-mapping variants; external-LD mismatch remains relevant. | Functional priors plus SuSiE PIPs/credible sets. | Genome-wide preprocessing adds substantial data, compute and dependency requirements; per-locus execution alone is not PolyFun. |

## Implemented p-LD\* definition

The literature treats LD pruning as a useful prioritization heuristic rather than a Bayesian fine-mapping model; the implemented project-specific algorithm is:

1. For each region and each configured threshold \(t\), independently sort variants by increasing P-value, with deterministic tie-breaking by chromosome, position and allele-aware variant key.
2. Visit variants in that order. Select a variant if its \(r^2\) with **every already selected variant** is strictly less than \(t\); otherwise exclude it and record a selected variant that caused exclusion.
3. Repeat from the full ranked list for every threshold, producing separate results. Thresholds always mean \(r^2\), not \(r\).

This is a greedy representative-selection algorithm; it is not an estimate of causal probability and is order-dependent by design. The project owner selected `r^2` thresholds 0.1, 0.01, and 0.001; each threshold produces an independent result. The current p-LD implementation processes all variants in the region, while top-N is an option for P-value ranking. Recommended safe handling for a required missing/out-of-band comparison is to mark the region as incomplete/error and re-extract or redefine the region; do not silently treat it as zero LD.

## Pan-UKB and input implications

The project specification says GWASprocessing provides GRCh38 GWAS-VCF, with ALT as effect allele and fields including `ES`, `SE`, `LP`, `AF`, `SS`, `EZ`, and provenance annotations. It also specifies Pan-UKB LD resources on GRCh37, sparse upper-triangular storage, covariate adjustment and a 10 Mb stored radius. These are treated as the intended input/reference contracts, not re-established by the fine-mapping papers reviewed here.

`ldcov` documentation describes extracting covariate-adjusted Pan-UKB LD, producing variant metadata and matrix outputs, and handling swapped alleles in its supported workflows ([official README](https://github.com/mkanai/ldcov/blob/92b8aad6b914ace1dad65eb119bae301f9448836/README.md)). It is an LD extraction utility, not a method that makes LD mismatch harmless. Verify its exact interface against the selected release and regional index before integrating it. The required GRCh38→GRCh37 mapping must be allele-aware; matching by rsID or coordinate alone cannot validate orientation. LD outside the stored band, absent variants, and ambiguous mappings must stay distinguishable from measured near-zero LD.

`SS` may vary across variants, while summary-statistic models can require a study or locus sample-size input. The implemented wrapper preserves per-variant values and uses the user-selected minimum non-missing regional `SS` unless an explicit effective-N override is supplied. This scalar remains an approximation, particularly for binary traits and meta-analyses, and requires representative-data validation. Likewise, a single ancestry-specific LD panel must not be inferred automatically from a GWAS ancestry label; the selected population, rationale, reference release and reference sample size belong in run provenance.

## Architecture and implementation recommendations

The literature supports the requested separation of concerns:

- **R analysis package** for GWAS-VCF validation, QC, allele-aware variant matching, loci/regions, LD validation/alignment, prioritization, SuSiE methods and result schemas. These functions should also be callable outside a workflow.
- **Nextflow** for task orchestration at approximately GWAS × region × method granularity, including SLURM resources, logs, provenance and aggregation. Share one validated/cached LD extraction across RSS and rRSS.
- **Python only for supporting components** whose established tools require it, such as a validated `ldcov` integration and later PolyFun; do not move the statistical core solely for orchestration.
- **Versioned configuration/input contracts** and preflight checks before expensive jobs; pin R/Python dependencies and record all data/reference versions.
- **Phase order**: P-value prioritization → p-LD\* → SuSiE-RSS → SuSiE-rRSS. Validate each before starting CARMA or PolyFun. Do not add FINEMAP, ABF or RSparsePro to implementation scope.

This architecture is a recommendation consistent with the user requirements, not a claim that workflow choice changes statistical performance.

## Confirmed decisions and remaining evidence gaps

Confirmed and implemented:

1. **p-LD\* defaults:** \(r^2 = 0.1, 0.01, 0.001\), analyzed independently; top-N applies only to P-value prioritization.
2. **`susieR` pin:** tag 2.0.0 at commit `4134b03dceb0e7c0e1a8885a926781ba4b0d6743`; synthetic tests exercise standard RSS and robust RSS diagnostics.
3. **Sample-size approximation:** minimum non-missing regional `SS`, with an explicit effective-N override. This is not an exact likelihood for heterogeneous sample sizes; binary-trait `SS`/`NC` interpretation needs validation.

Remaining evidence gaps:

1. **Locus/gene region rules:** genome-wide significance, merging distance, extension, gene annotation/build and cis-window need evidence-backed defaults and per-run overrides. Cis/trans labeling must remain separate from fine-mapping.
2. **Liftover and ambiguous alleles:** select a chain/map source and define indel, strand, palindromic, multiallelic and ambiguous handling; retain every mapping/exclusion reason.
3. **LD QC and repair:** define configured checks for symmetry, diagonal, missingness, correlation range, PSD and any explicit regularization/repair. Never repair silently.
4. **Population and multi-ancestry policy:** choose LD population explicitly, record the rationale, and do not treat a multi-ancestry meta-analysis as if one panel were automatically matched.
5. **Phase 2 evidence and validation:** obtain the CARMA paper text and review inaccessible supplements before selecting or implementing CARMA/PolyFun defaults. Benchmark actual representative regions before promising resource requirements.

## Additional sources reviewed

- Zou Y, Carbonetto P, Wang G, Stephens M. 2022. “Fine-mapping from summary data with the ‘Sum of Single Effects’ model.” *PLOS Genetics* 18:e1010299. [DOI](https://doi.org/10.1371/journal.pgen.1010299); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9337707/).
- Weissbrod O et al. 2020. “Functionally informed fine-mapping and polygenic localization of complex trait heritability.” *Nature Genetics* 52:1355–1363. [DOI](https://doi.org/10.1038/s41588-020-00735-5); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7710571/).
- Yang Z et al. 2023. “CARMA is a new Bayesian model for fine-mapping in genome-wide association meta-analyses.” *Nature Genetics* 55:1057–1065. [DOI](https://doi.org/10.1038/s41588-023-01392-0); [official package documentation](https://github.com/ZikunY/CARMA).
- Zhang W et al. 2026. “Robust fine-mapping in the presence of linkage disequilibrium mismatch.” *Nature Genetics* 58:2993–3003. [DOI](https://doi.org/10.1038/s41588-026-02780-y); [RSparsePro code](https://github.com/zhwm/RSparsePro_LD).
- `susieR` version-specific source: [tag 2.0.0 `susie_rss` documentation](https://github.com/stephenslab/susieR/blob/4134b03dceb0e7c0e1a8885a926781ba4b0d6743/R/susie.R), [RSS mismatch vignette](https://github.com/stephenslab/susieR/blob/4134b03dceb0e7c0e1a8885a926781ba4b0d6743/vignettes/rss_mismatch.Rmd), [CRAN package history](https://cran.r-project.org/package=susieR).
- GeneHackman implementation reference: [MRCIEU/GeneHackman `R/finemap.R`](https://github.com/MRCIEU/GeneHackman/blob/master/R/finemap.R). It illustrates ordinary SuSiE-RSS and ancestry-specific LD integration, not the newer robust settings.
- Pan-UKB extraction reference: [`ldcov` README](https://github.com/mkanai/ldcov/blob/92b8aad6b914ace1dad65eb119bae301f9448836/README.md).

