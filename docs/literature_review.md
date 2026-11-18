# Fine-mapping literature review

**Status:** Literature findings are complete for the locally supplied full texts and accessible primary sources. Several linked supplementary files and the CARMA article full text could not be retrieved; those gaps are recorded below. Recommendations are provisional where they depend on those gaps or on project choices.

## Scope and evidence

The review covers the five papers supplied in `docs/papers/`, the primary SuSiE-RSS and PolyFun papers, current `susieR` documentation, CARMA implementation documentation, RSparsePro as a relevant LD-mismatch comparator, and GeneHackman and `ldcov` as implementation references. Literature describes statistical behavior; the Pan-UKB release facts and GWAS-VCF field contract below are project requirements and still need release-specific preflight verification.

### Supplied-paper inventory

| File | Paper |
|---|---|
| [2021.09.03.21262975v1.full.pdf](./papers/2021.09.03.21262975v1.full.pdf) | Kanai et al. (2021), “Insights from complex trait fine-mapping across diverse populations,” *medRxiv*, DOI [10.1101/2021.09.03.21262975](https://doi.org/10.1101/2021.09.03.21262975). |
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

Yang et al. (2023), “CARMA is a new Bayesian model for fine-mapping in genome-wide association meta-analyses,” *Nature Genetics* 55:1057–1065, DOI [10.1038/s41588-023-01392-0](https://doi.org/10.1038/s41588-023-01392-0), is the primary paper. The publisher version-of-record full text was not accessible via open-access endpoints; implementation and parameter semantics were verified directly against the pinned source repository (`ZikunY/CARMA@856f0fb48cc42f1ad895600f0a4d80ad8451c2ff`), [official CARMA manual](https://github.com/ZikunY/CARMA/blob/master/man/CARMA.Rd), `CARMA_demo.pdf`, and underlying R/C++ routines (`ind_Normal_fixed_sigma_marginal`, `outlier_ind_Normal_marginal`). CARMA operates on summary Z-statistics and signed LD correlations (`rho.index = 0.99`, `EM.dist = "Logistic"`, `Max.Model.Dim = 2e5`, `all.iter = 3`, `all.inner.iter = 10`, `outlier.BF.index = 1 / 3.2`), producing variant-level PIPs, purity-filtered credible sets, and discrepancy-outlier indices. Do not infer that outlier removal resolves all forms of population, allele, or reference-size mismatch.

### PolyFun + SuSiE

Weissbrod et al. (2020), “Functionally informed fine-mapping and polygenic localization of complex trait heritability,” *Nature Genetics* 52:1355–1363, DOI [10.1038/s41588-020-00735-5](https://doi.org/10.1038/s41588-020-00735-5) (open-access full text reviewed via Europe PMC [PMC7710571](https://europepmc.org/articles/PMC7710571)), introduces PolyFun in two operational modes: (1) trait-specific genome-wide regularized S-LDSC + L2-binned cross-fitting, and (2) **Approach 1 (parametric / precomputed meta-analyzed UK Biobank functional priors, `snpvar_meta`)**, which requires no trait-specific S-LDSC refitting (`polyfun.py` / `finemapper.py` at pinned commit `omerwe/polyfun@8b17b4ec87fde11f8d216c049c08d1420cd25de6`). Precomputed autosomal GRCh37 `snpvar_meta` priors across ~19.48M common/low-frequency UKB variants (187 Baseline-LF2.2.UKB annotations meta-analyzed across 15 UKB traits, with $q = 100$ truncation floor $\max_j \sigma_j^2 / 100$) are matched allele-directionally on lifted GRCh37 `(CHR, BP, REF, ALT)` coordinates and passed into `susieR::susie_rss(..., prior_weights = snpvar / sum(snpvar))`. Unmatched (e.g. rare) variants receive the locus baseline floor (`missing_prior_policy = "baseline_floor"`) and are explicitly counted in `qc.json`.

## Method comparison

Standalone `p-value` output generation was retired because Stage 1 (`scripts/screen_gwas.py`) already screens and identifies genome-wide significant leads ($P \le 5 \times 10^{-8}$). All five methods below run by default (`--methods "p-ld susie-rss susie-rrss carma polyfun"`), or any user-specified subset can be selected.

| Method | Model / inputs | LD assumptions and mismatch | Outputs | Compute profile and principal limitations |
|---|---|---|---|---|
| **p-LD\*** (`p-ld`) | Greedy P-value-ordered pruning using pairwise \(r^2\); input statistics plus regional LD. | Requires usable, correctly matched \(r^2\) for the comparisons made. Missing/out-of-band LD is unknown, not zero. Sensitive to threshold and variant set. | Selected variants per \(r^2\) threshold (`0.1`, `0.01`, `0.001`), rank and pruning/exclusion reason. **No PIPs or credible sets.** | Fast (~4–14s per 3 Mb locus). Heuristic; thresholds are not posterior uncertainty. |
| **SuSiE-RSS** (`susie-rss`) | Bayesian sum-of-single-effects regression; aligned Z, LD matrix, sample size, `L = 10`, `max_iter = 50`. | External LD approximates in-sample LD; requires compatible variants, allele coding, ancestry and sufficient accuracy. Mismatch can distort PIPs and coverage. | PIPs, per-effect credible sets, purity and model/convergence diagnostics, plus `coloc.susie`-ready `alpha`/`lbf_variable`. | ~35–50s per 3 Mb locus with in-place C `XtX` scaling (`<6 GB` R heap). Sensitive to severe multi-ancestry LD mismatch. |
| **SuSiE-rRSS** (`susie-rrss`) | SuSiE-RSS with finite-reference (`R_finite`) and LD-mismatch (`R_mismatch = "eb_mix"`) modeling in pinned `susieR` 2.0.0. | Uses the same aligned Z/R requirements; `R_finite` and `R_mismatch` address finite-panel and regional/variant mismatch uncertainty. | SuSiE outputs plus robust settings, mismatch/reliability diagnostics and `SER_fallback`/`R_reliability_flag`. | ~55–185s per 3 Mb locus. Resolves non-convergence on multi-ancestry `ALL -> EUR` loci while flagging severe mismatch. |
| **CARMA** (`carma`) | Shotgun stochastic search spike-and-slab model (`ZikunY/CARMA`); aligned Z and signed LD `R`, Poisson model-size prior, `outlier.switch = TRUE`. | Explicit Bayes-factor (`BF < 1/3.2`) LD-discrepancy outlier detection on candidate causal configurations. | Variant PIPs, 95% credible sets (`min_abs_corr >= 0.5`), and flagged outlier variant IDs/indices. | ~20–28s per 3 Mb locus (`<4.3 GB` R heap) via vectorised LAPACK batch marginal-likelihood evaluation and sparse prior handling. |
| **PolyFun + SuSiE** (`polyfun`) | Precomputed UKB meta-analyzed `snpvar_meta` functional priors (Weissbrod et al. 2020 Approach 1) + `susieR::susie_rss(..., prior_weights = ...)`. | Requires GRCh37 `(CHR, BP, REF, ALT)` match to UKB Baseline-LF2.2 annotations (autosomes 1–22); external-LD mismatch caveats match SuSiE-RSS. | Prior-weighted PIPs, 95% credible sets, per-SNP `SNPVAR` & `PRIOR_MATCHED`, and `coloc.susie`-ready matrices. | ~55–218s per 3 Mb locus (`<5.9 GB` R heap). Concentrates posterior mass onto functionally enriched coding/regulatory variants within high-LD blocks. |

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

## Locus clustering, window defaults, significance gates, and MHC policy (resolved 17 Nov 2026)

Full-text review of all six supplied PDFs alongside open-access full texts retrieved via Europe PMC — Sun et al. (2023, UKB-PPP, *Nature*, [PMC10567551](https://europepmc.org/article/PMC/PMC10567551)), Eldjarn et al. (2023, *Nature*, [PMC10567571](https://europepmc.org/articles/PMC10567571)), Wallace (2021, `coloc.susie`, *PLOS Genetics*, [PMC8504726](https://europepmc.org/articles/PMC8504726)), Zou et al. (2022, `SuSiE-RSS`, [PMC9337707](https://europepmc.org/articles/PMC9337707)), Benner et al. (2016, `FINEMAP`, [PMC4866522](https://europepmc.org/articles/PMC4866522)), Hutchinson, Watson & Wallace (2020, [PMC7179948](https://europepmc.org/articles/PMC7179948)), Weissbrod et al. (2020, `PolyFun`, [PMC7710571](https://europepmc.org/articles/PMC7710571)), Kerimov et al. (2021, `eQTL Catalogue`, [PMC8423625](https://europepmc.org/articles/PMC8423625)), Barbeira et al. (2021, GTEx v8, [PMC7836161](https://europepmc.org/articles/PMC7836161)), and Broekema et al. (2020, [PMC7014684](https://europepmc.org/articles/PMC7014684)) — plus source-code audits of `mkanai/finemapping-pipeline` and `MRCIEU/GeneHackman` establishes four production rules:

1. **Regional significance gate ($P \le 5 \times 10^{-8}$, `LP >= 7.30103` by default across all region types):**
   - Every reviewed GWAS/pQTL method, benchmark, and biobank atlas (Schaid et al. 2018 p. 2; Benner et al. 2016 $p_0 = 0$; Kanai et al. 2021; Sun et al. 2023; Eldjarn et al. 2023; Li & Zhou 2025 pp. 2–3; Chen et al. 2026) requires genome-wide significance before fine-mapping both `cis` and `trans` loci.
   - In pinned `susieR` 2.0.0, `check_null_threshold = 0` and `susie_get_pip(..., prune_by_cs = FALSE)` allow sub-threshold noise or minor external-LD discrepancies to emit $\sum \text{PIP} \ge 1$ or false-positive credible sets in weak/null windows (Zou et al. 2022; Hutchinson et al. 2020). Wallace (2021, Fig. 2) shows `coloc.susie` performs 0 comparisons when either trait has 0 credible sets and recommends falling back to single-variant `coloc.abf` from marginal `ES`/`SE` rather than running `SuSiE` on underpowered loci.
   - Secondary credible sets inside an activated ($P \le 5 \times 10^{-8}$) locus must **never** be filtered by marginal P-value so opposite-LD masked signals are retained (Benner et al. 2016; Kanai et al. 2021 p. 9; Li & Zhou 2025 p. 4). An explicit provenance-recorded override (`--cis-significance-lp`) supports specialist cis-QTL designs.
2. **Greedy 3 Mb (`±1.5 Mb`) lead-centred windowing + overlap merge (vs static LD blocks or `GeneHackman`):**
   - Kanai et al. (2021, pp. 3, 9, 19; `make_finemap_inputs.py`) rank $P < 5 \times 10^{-8}$ variants by $\chi^2$ descending, place a $\pm 1.5\text{ Mb}$ (3 Mb) window on the top variant, drop other significant variants inside that window from lead selection, repeat, and merge overlapping intervals — capturing distal non-coding enhancers 0.6–1.2 Mb away (`PHACTR1`/`EDN1`, `GATA3`, `8q24`) and long-range negative-LD masking (`FGFR4`, `CETP`).
   - Static haplotype/LD-block boundaries are rejected per Schaid et al. (2018, p. 7) because arbitrary block edges truncate correlated variants near boundaries. `MRCIEU/GeneHackman`'s `plink --clump` + `±1 Mb` approach is rejected because it relies on 1000G rsIDs, fails to merge overlapping windows in single-ancestry mode (`R/finemap.R`), and lacks MHC handling.
   - For protein GWAS, target cis regions remain anchored to `gene start − 1 Mb` to `gene end + 1 Mb` (`region_type = "cis_gene"`), harmonising boundaries across Olink, SomaScan, and GTEx (Sun et al. 2023). A 9.5 Mb valley-split ceiling prevents merged chains from exceeding Pan-UKB's 10 Mb stored LD radius.
3. **Single-region MHC (`chr6:25,500,000–34,000,000` GRCh38) + `mhc_region` flags + `--exclude-mhc`:**
   - Official NCBI GRC `assembly_regions.txt` defines the core MHC as `chr6:28,477,797–33,448,354` (GRCh37.p13) and `chr6:28,510,120–33,480,577` (GRCh38.p14). Following Sun et al. (2023), treating extended HLA (`chr6:25,500,000–34,000,000` GRCh38, 8.50 Mb) as **1 single locus** captures extended `HIST1` + HLA class I/II/III LD while staying inside Pan-UKB's 10 Mb BlockMatrix radius (whereas Kanai's prose `25–36 Mb` = 11 Mb exceeds the 10 Mb band). `MRCIEU/GeneHackman` has no MHC handling and slices HLA into overlapping 2 Mb clumps.
4. **Unobtained full texts (explicitly excluded from citation):**
   - Full texts were not accessible in `docs/papers/` or Europe PMC Open Access (`OPEN_ACCESS:y`) for: Wang et al. (2020, *JRSS-B* `10.1111/rssb.12388`), Fernandez-Rozadilla et al. (2023, *Nature Genetics* `10.1038/s41588-022-01222-9`), Berisa & Pickrell (2016, *Bioinformatics* `10.1093/bioinformatics/btv546`), Ferkingstad et al. (2021, `10.1038/s41588-021-00978-w`), Pietzner et al. (2021, `10.1126/science.abj1541`), GTEx Consortium (2020, *Science* `10.1126/science.aaz1776`), Hormozdiari et al. (2014, `CAVIAR`), Yang et al. (2023, `CARMA`), Zhang et al. (2026, `RSparsePro` VoR), Cui et al. (2024, `SuSiE-inf`), Yuan et al. (2024, `SuSiEx`), Gao & Zhou (2024, `MESuSiE`), and publisher supplementary files for Kanai (2021), Sun (2023), Benner (2017), and Schaid (2018).

## Confirmed decisions and remaining evidence gaps

Confirmed:

1. **p-LD\* defaults:** \(r^2 = 0.1, 0.01, 0.001\), analyzed independently; top-N applies only to P-value prioritization.
2. **`susieR` pin:** tag 2.0.0 at commit `4134b03dceb0e7c0e1a8885a926781ba4b0d6743`; synthetic tests exercise standard RSS and robust RSS diagnostics.
3. **Sample-size approximation:** minimum non-missing regional `SS`, with an explicit effective-N override. This is not an exact likelihood for heterogeneous sample sizes; binary-trait `SS`/`NC` interpretation needs validation.
4. **Significance gate, windows, and MHC:** default $P \le 5 \times 10^{-8}$ gate across all regions (opt-in `--cis-significance-lp`), Kanai greedy `±1.5 Mb` (3 Mb) merged loci (9.5 Mb ceiling), gene-anchored `cis_gene` (`±1 Mb`), and single merged extended MHC (`chr6:25.5–34.0 Mb` GRCh38) with `mhc_region` / `mhc_core_grc` flags and `--exclude-mhc`.

Remaining evidence gaps:

1. **Liftover and ambiguous alleles:** retain every mapping/exclusion reason (`no_chain_mapping`, `ref_mismatch_after_liftover`, `no_exact_panel_allele_match`, `duplicate_panel_variant`) in `qc.json` and monitor attrition across studies.
2. **LD QC and repair:** validate symmetry, diagonal, finite values, and correlation range on every locus (skipping $O(M^3)$ full `eigen()` when $M > 5,000$ such as single-region MHC). Never repair silently.
3. **Population and multi-ancestry policy:** choose LD population explicitly among the 6 Pan-UKB panels (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`), record its sample size and rationale, and do not treat a multi-ancestry meta-analysis as if one panel were automatically matched.
4. **CARMA and PolyFun integration:** `CARMA` (`ZikunY/CARMA@856f0fb4`) and `PolyFun + SuSiE` (`omerwe/polyfun@8b17b4ec` UKB `snpvar_meta` priors) are integrated alongside `p-ld`, `susie-rss`, and `susie-rrss` as default pipeline methods and validated across Locke 2015 BMI and Sun 2023 `NPPB`. Note that Yang et al. (2023) publisher VoR text remained paywalled during review; CARMA parameter defaults follow upstream `ZikunY/CARMA` package documentation and source code.

## Additional sources reviewed

- Zou Y, Carbonetto P, Wang G, Stephens M. 2022. “Fine-mapping from summary data with the ‘Sum of Single Effects’ model.” *PLOS Genetics* 18:e1010299. [DOI](https://doi.org/10.1371/journal.pgen.1010299); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9337707/).
- Weissbrod O et al. 2020. “Functionally informed fine-mapping and polygenic localization of complex trait heritability.” *Nature Genetics* 52:1355–1363. [DOI](https://doi.org/10.1038/s41588-020-00735-5); [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7710571/).
- Yang Z et al. 2023. “CARMA is a new Bayesian model for fine-mapping in genome-wide association meta-analyses.” *Nature Genetics* 55:1057–1065. [DOI](https://doi.org/10.1038/s41588-023-01392-0); [official package documentation](https://github.com/ZikunY/CARMA).
- Zhang W et al. 2026. “Robust fine-mapping in the presence of linkage disequilibrium mismatch.” *Nature Genetics* 58:2993–3003. [DOI](https://doi.org/10.1038/s41588-026-02780-y); [RSparsePro code](https://github.com/zhwm/RSparsePro_LD).
- `susieR` version-specific source: [tag 2.0.0 `susie_rss` documentation](https://github.com/stephenslab/susieR/blob/4134b03dceb0e7c0e1a8885a926781ba4b0d6743/R/susie.R), [RSS mismatch vignette](https://github.com/stephenslab/susieR/blob/4134b03dceb0e7c0e1a8885a926781ba4b0d6743/vignettes/rss_mismatch.Rmd), [CRAN package history](https://cran.r-project.org/package=susieR).
- GeneHackman implementation reference: [MRCIEU/GeneHackman `R/finemap.R`](https://github.com/MRCIEU/GeneHackman/blob/master/R/finemap.R). It illustrates ordinary SuSiE-RSS and ancestry-specific LD integration, not the newer robust settings.
- Pan-UKB extraction reference: [`ldcov` README](https://github.com/mkanai/ldcov/blob/92b8aad6b914ace1dad65eb119bae301f9448836/README.md).

