# Architecture and data contracts

**Status:** Phase A–C design, methods contracts, and local production implementation finalised (17 Nov 2026); targets >10,000-GWAS storage efficiency, evidence-backed region/MHC defaults, idempotent batch execution over `GWASprocessing` outputs, self-contained HTML reporting, and container/SLURM portability. Reusable epidemiological and statistical fine-mapping defaults live in `~/.claude/skills/methods/references/fine-mapping.md`; this file specifies the `GWASfinemapping` software architecture and input/output contracts. See [the Pan-UKB LD guide](panukb-ld.md) and [the reference-data guide](reference-data.md).

## System boundaries

`GWASprocessing` is the upstream producer of standardised GRCh38 GWAS-VCF files. `GWASfinemapping` consumes those files; it does not make `GWASprocessing` a runtime dependency or repeat its harmonisation. This repository owns its input validation, allele-aware liftover and Pan-UKB reference matching, statistical methods, Nextflow orchestration, storage-light output contract, container environment, test suite, and benchmark.

The statistical core lives in the `gwasfinemapping` R package (`R/`) so functions can be invoked interactively, in unit tests, and from workflow runners. Nextflow (`finemapping.nf`) and `scripts/finemap.sh` orchestrate local and SLURM execution across single GWAS or multi-GWAS manifests. Python (`scripts/run_finemapping.py`, `scripts/reference_data.py`, `scripts/report_finemapping.py`) handles indexed VCF extraction (`bcftools`), PyArrow regional slicing of installed Pan-UKB variant indexes, `ldcov` S3 BlockMatrix extraction, scratch lifecycle cleanup, and self-contained HTML report rendering.

## Components

| Component | Responsibility |
|---|---|
| `scripts/source.sh` + `scripts/lib/env.sh` | Bootstrap micromamba + Python tools (`Hail`, `ldcov[s3]`) + isolated R library (`gwasfinemapping`, `susieR` 2.0.0) + Apptainer `.sif` image; redirect work/cache/reference paths out of cloud-synced folders. |
| `scripts/install_reference_data.sh` | Install and verify UCSC `hg38ToHg19.over.chain.gz`, Ensembl GRCh37 FASTA + `.fai`, UniProt GRCh38 gene coordinates (`usefuldata`), and all 6 whole-population Pan-UKB GRCh37 Parquet variant indexes (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`). |
| `scripts/prepare_manifest.py` + `finemapping.nf` | Discover standardised GWAS across `GWASprocessing` study trees (`--input-dir`), explicit manifests (`--manifest`), or single VCFs (`--vcf`); support `--only-study`, `--row-limit`, `--dry-run`, `--force`, and `--report true\|false`; stage safe per-row JSON configs; track completed runs idempotently in `<outdir>/.finemap_state.tsv`; and isolate per-row failures with `<outdir>/<study>/logs/` logs. |
| Region builder (`scripts/run_finemapping.py`) | Enforce the default $P \le 5 \times 10^{-8}$ (`LP >= 7.30103`) significance gate across all region types; construct target `cis_gene` (`gene start − 1 Mb` to `gene end + 1 Mb`), greedy lead-centred 3 Mb (`±1.5 Mb`) merged genome-wide loci (with a 9.5 Mb valley-split ceiling for Pan-UKB's 10 Mb LD band), and a single merged extended MHC locus (`chr6:25,500,000–34,000,000` GRCh38, or skip when `--exclude-mhc` is set). |
| Scratch LD & liftover stage | Slice the installed population Parquet index via PyArrow (`cut_regional_index`), lift GRCh38 variants to GRCh37 with swapped-allele indel pre-normalisation and `bcftools norm`, extract Pan-UKB LD via `ldcov`, pass aligned float64 LD to R via binary `.npy` (no ASCII `matrix.tsv.gz`), and delete per-locus scratch files immediately on locus completion unless `--keep-intermediates` is requested. |
| R statistical package (`R/`) | Validate GWAS records (`validate_gwas_vcf_records`) and LD matrices (`validate_ld_matrix`, skipping $O(M^3)$ dense `eigen()` when $M > 5,000$); run deterministic `p-value` and `p-ld` ($r^2 < 0.1, 0.01, 0.001$) prioritisation plus Bayesian `SuSiE-RSS` and `SuSiE-rRSS` (`susieR` 2.0.0) with `GWAS_VARIANT_ID`-named `alpha`, `lbf_variable`, and `pip` objects ready for `coloc::coloc.susie()` and `coloc::coloc.bf_bf()`. |
| HTML reporting (`render_report.R` + `render_manifest_report.py`) | Render one self-contained `<study_id>/<gwas_id>/report.html` per GWAS directly from saved `.rda` + `.json` outputs (also callable on demand after `--report false` batch runs) plus a batch-level `manifest.html`. |
| Validation & benchmark (`scripts/test.sh` + `scripts/benchmark.sh`) | Run the 4-stage offline test suite (`pytest`, `testthat`, `R CMD check`, Nextflow integration) and the reproducible real-data benchmark on Locke 2015 BMI (`GWAScatalog_GCST002783`) and Sun 2023 UKB-PPP `NPPB` (`P16860_OID20049`) against `benchmark/expected/*.json` (plus `--offline` and `--server-panukb` modes). |

## GWAS-VCF input contract

Inputs are bgzipped, tabix-indexed GRCh38 GWAS-VCF files (`*.vcf.gz` + `.tbi`) produced by `GWASprocessing`, with `ALT` as the effect allele:

| Field | Contract |
|---|---|
| `CHROM`, `POS`, `REF`, `ALT`, `ID` | Required GRCh38 coordinates and alleles (`1`–`22`, `X`, `Y`); canonical `GWAS_VARIANT_ID` is `chr:pos:ref:alt` on GRCh38. |
| `FORMAT/ES`, `FORMAT/SE` | Effect size relative to `ALT` and positive standard error; finite numeric values required. |
| `FORMAT/LP` | $-\log_{10}(P)$; finite non-negative values required (`LP >= 7.30103` corresponds to $P \le 5 \times 10^{-8}$). |
| `FORMAT/SS` | Per-variant sample size; regional `SuSiE` uses minimum non-missing regional `SS` unless `--effective-n` overrides it. |
| `FORMAT/AF`, `FORMAT/AFS`, `INFO/HM` | Extracted and preserved in per-locus variant tables when present (`NA` if absent in source VCF) for downstream Mendelian randomisation palindromic/MAF checks and `coloc`. |

## Region construction and MHC contract

1. **Significance gate (`--significance-lp`, default `7.30103`, $P \le 5 \times 10^{-8}$):**
   - Applied uniformly by default to genome-wide loci, `cis_gene` windows, and the MHC before LD extraction or method execution.
   - Optional explicit override `--cis-significance-lp` (recorded in `provenance.json`) allows specialist cis-QTL thresholds when pre-specified. If `cis_gene` has no variant meeting its gate, LD extraction and fine-mapping are skipped and `qc.json` records `status = "skipped_below_significance_threshold"` alongside `max_lp`.
   - Secondary credible sets inside an activated locus are **never** filtered by marginal P-value.
2. **Target `cis_gene` window (`--gene-window-bp`, default `1,000,000`):**
   - When `--uniprot-id` or `--gene-id` resolves uniquely in `gene_regions_GRCh38.tsv`, constructs `[max(1, start - 1e6), end + 1e6]` (`region_type = "cis_gene"`).
3. **Genome-wide lead-centred windows (`--region-flank-bp`, default `1,500,000` = 3 Mb window):**
   - Greedy positional lead selection following Kanai et al. (2021): rank significant non-MHC variants by `LP` descending, centre `[lead - 1.5 Mb, lead + 1.5 Mb]` on the top variant, mask all other significant variants inside that 3 Mb window from further lead selection, repeat until exhausted, and **merge overlapping intervals**.
   - **Pan-UKB 10 Mb radius guard (`--max-region-width-bp`, default `9,500,000`):** If chained merges exceed 9.5 Mb on GRCh38/GRCh37, split at the deepest $-\log_{10} P$ valley between leads and record `split_at_panukb_radius = true` in `qc.json`.
4. **MHC region (`chr6:25,500,000–34,000,000` GRCh38; official GRC core `chr6:28,510,120–33,480,577` GRCh38 / `chr6:28,477,797–33,448,354` GRCh37):**
   - Default behaviour: collapse all significant MHC variants into **1 single `chr6:25,500,000–34,000,000` region** (`region_type = "mhc"`, 8.50 Mb, strictly within Pan-UKB's 10 Mb stored LD radius), truncate adjacent non-MHC flanks at `25,499,999` / `34,000,001`, and set `mhc_region = TRUE` and `mhc_core_grc` on locus summaries, variant tables, and credible sets.
   - Pipeline option `--exclude-mhc` (`params.exclude_mhc`): records discovered MHC hits in `qc.json` (`status = "excluded_by_option"`) and skips MHC LD extraction and fine-mapping.

## Storage-light output contract (>10,000 GWAS scale)

To keep 10,000 fine-mapped GWAS within ~200–240 GB (~20–24 MB per GWAS instead of ~10 GB per GWAS of dense LD intermediates), each published `<output_root>/<study_id>/<gwas_id>/` directory retains **only**:

```text
<study_id>/<gwas_id>/
├── metadata.json          # Study/trait/gene/panel summary & status
├── provenance.json        # Software/git versions, CLI parameters, reference SHA-256s, Pan-UKB sample size
├── qc.json                # Complete locus manifest, liftover/panel attrition counts, LD validation, CS & mismatch flags
├── report.html            # Self-contained QC & locus visualisation report (base64 plots, zero external assets)
└── methods/
    ├── p-value.rda        # Standalone xz-compressed R object `p_value`
    ├── p-ld.rda           # Standalone xz-compressed R object `p_ld` (r2 = 0.1, 0.01, 0.001)
    ├── susie-rss.rda      # Standalone xz-compressed R object `susie_rss` (coloc.susie-ready `fit`)
    └── susie-rrss.rda     # Standalone xz-compressed R object `susie_rrss` (R_finite + R_mismatch = "eb_mix")
```

All per-locus intermediate directories (`discovery/*.tsv`, `gene_region/*.tsv`, `loci/<id>/inputs/`, `loci/<id>/mapping/`, `loci/<id>/ld/*.{npz,npy,tsv}`, `loci/<id>/methods/`) are written to scratch space and deleted automatically after aggregation and report rendering (unless `--keep-intermediates` is explicitly passed for debugging).

### `coloc.susie` & `coloc.bf_bf` compatibility contract

In both `susie-rss.rda` and `susie-rrss.rda`, every region's `fit` object (`susie_rss$regions[[region_id]]$result$fit`) sets:
- `colnames(fit$alpha)`, `colnames(fit$lbf_variable)`, `colnames(fit$mu)`, `colnames(fit$mu2)`, and `names(fit$pip)` to canonical GRCh38 `GWAS_VARIANT_ID` (`chr:pos:ref:alt`);
- `fit$sets$cs`, `fit$sets$cs_index`, and `fit$sets$purity` from `susieR` (`min_abs_corr = 0.5`, `coverage = 0.95`);
- Per-variant `ID` (rsID), `PANEL_VARIANT_ID` (GRCh37), `ES`, `SE`, `LP`, `SS`, `AF`, `AFS`, `HM`, `mhc_region`, and `mhc_core_grc` in `result$variants`.

Downstream `coloc::coloc.susie(fit1, fit2)` and `coloc::coloc.bf_bf()` therefore run directly on the saved `.rda` files with no LD matrix required.

## Confirmed decisions

1. **Standalone per-method `.rda` layout (`Option 1B`):** Four self-contained `xz`-compressed `.rda` files per GWAS (`p-value.rda`, `p-ld.rda`, `susie-rss.rda`, `susie-rrss.rda`), working identically for single complex-trait GWAS (such as BMI or colon and rectal cancer) and multi-thousand-protein studies (such as UKB-PPP or deCODE).
2. **Significance gate:** Default `LP >= 7.30103` ($P \le 5 \times 10^{-8}$) across all region modes, with explicit opt-in `--cis-significance-lp`.
3. **Window rules:** Kanai greedy `±1.5 Mb` (3 Mb) merged loci (9.5 Mb valley-split ceiling) + gene-anchored `cis_gene` (`±1 Mb`) + single merged extended MHC (`chr6:25.5–34.0 Mb` GRCh38) with `mhc_region` flags and `--exclude-mhc`.
4. **All 6 Pan-UKB LD populations supported automatically:** `AFR` (`6,637`), `AMR` (`982`), `CSA` (`8,876`), `EAS` (`2,709`), `EUR` (`420,542`), and `MID` (`1,599`), verified from the Pan-UKB Hail table globals.
5. **`susieR` pin:** GitHub tag `2.0.0` (`4134b03dceb0e7c0e1a8885a926781ba4b0d6743`).
