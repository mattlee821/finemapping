# GWASfinemapping

Scalable statistical fine-mapping and variant prioritisation for standardised GRCh38 GWAS-VCF summary statistics (designed for both individual complex-trait GWAS and >10,000-trait proteogenomic studies) against the six public Pan-UKB LD panels (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`).

`GWASfinemapping` consumes the standardised GWAS-VCF outputs of `GWASprocessing` without runtime coupling: see [docs/architecture.md](docs/architecture.md), [docs/literature_review.md](docs/literature_review.md), [docs/panukb-ld.md](docs/panukb-ld.md), and [docs/reference-data.md](docs/reference-data.md).

## Key features

- **Cascading 4-stage architecture with upfront screen & shared LD prep:**
  1. `scripts/screen_gwas.py`: fast $P \le 5 \times 10^{-8}$ (`LP >= 7.30103`) significance screen across all region types (early-exiting GWAS with 0 genome-wide significant variants in ~30s), building Kanai et al. (2021) greedy `±1.5 Mb` (**3 Mb**) merged lead windows (9.5 Mb ceiling), gene-anchored `cis_gene` (`gene ± 1 Mb`) windows, and a single merged extended MHC locus (`chr6:25.5–34.0 Mb` GRCh38, or skipped via `--exclude-mhc`).
  2. `scripts/prep_locus.py` (and interactive `scripts/interactive_prep_locus.R`): shared per-locus GRCh38→GRCh37 liftover, FASTA-verified normalisation, Pan-UKB `VariantIndex` matching, and Fortran-order `.npy` LD extraction reused across all fine-mapping methods.
  3. Standalone, interactive method scripts (`scripts/methods/run_{p_ld,susie_rss,susie_rrss,carma,polyfun}.R`): runnable both line-by-line in RStudio/Positron and in parallel across Nextflow tasks.
  4. `scripts/finalise_gwas.py`: aggregates per-locus outputs into GWAS-level `.rda` / `.rds` method files, writes `{metadata,provenance,qc}.json` (including per-step timing & memory telemetry) and self-contained `report.html`, and removes multi-gigabyte scratch LD matrices immediately.
- **Five fine-mapping & prioritisation methods enabled by default (selectable via `--methods`):**
  - `p-ld`: deterministic greedy $r^2$-pruning ($r^2 < 0.1, 0.01, 0.001$) via a Selected-Push Sieve (`34×` faster than row-scan clumping).
  - `susie-rss`: summary-statistic `SuSiE-RSS` (`susieR` 2.0.0, zero-copy in-place $X^\top X = (n - 1)R$ scaling).
  - `susie-rrss`: mismatch-aware `SuSiE-rRSS` (`R_finite`, `R_mismatch = "eb_mix"`).
  - `carma`: outlier-robust Bayesian Shotgun fine-mapping (`ZikunY/CARMA`, Spike-slab / Cauchy priors + discrepancy Bayes-factor outlier detection).
  - `polyfun`: functionally informed `PolyFun + SuSiE` (Weissbrod et al. 2020 UK Biobank `baselineLF2.2.UKB` per-SNP heritability priors with $q = 100$ truncation, or custom S-LDSC priors).
  - Pass any space- or comma-separated subset with `--methods`, e.g. `--methods "p-ld susie-rrss polyfun"`.
- **Storage-light, `coloc`-ready outputs (~25–35 MB/GWAS with `xz`):** Retains standalone `xz`-compressed `.rda` (or `.rds`) files (`methods/{p-ld,susie-rss,susie-rrss,carma,polyfun}.rda` with `GWAS_VARIANT_ID`-named `alpha` and `lbf_variable` matrices ready for `coloc::coloc.susie()` and `coloc::coloc.bf_bf()`), `{metadata,provenance,qc}.json`, and `report.html`.
- **Portable local, container & SLURM setup:** One-command bootstrap (`scripts/source.sh`), reference installer (`scripts/install_reference_data.sh`), validation test suite (`scripts/test.sh`), and reproducible live-S3 + synthetic offline benchmark (`scripts/benchmark.sh`) matching `GWASprocessing`.

## Quickstart

```bash
# 1. Install pinned conda/Python/R environment (including susieR 2.0.0 and CARMA)
bash scripts/source.sh
source scripts/source.sh

# 2. Install & verify reference data once (~3.8 GB under ~/.cache/gwasfinemapping/reference)
bash scripts/install_reference_data.sh install
bash scripts/install_reference_data.sh check

# 3. Validate installation (unit tests, R CMD check, Nextflow DAG + live-S3 benchmark)
bash scripts/test.sh
bash scripts/benchmark.sh            # Locke 2015 BMI + Sun 2023 UKB-PPP NPPB live from Pan-UKB S3
bash scripts/benchmark.sh --offline  # zero-network synthetic benchmark (all 5 methods)

# 4. Run all 5 methods (default) or a specific subset locally or on SLURM
bash scripts/finemap.sh --input-dir /path/to/GWASprocessing/data/GWAS --ld-pop EUR --outdir results/ --dry-run
bash scripts/finemap.sh --input-dir /path/to/GWASprocessing/data/GWAS --methods "p-ld susie-rrss polyfun" --ld-pop EUR --outdir results/
sbatch scripts/finemap.sh --input-dir /path/to/GWASprocessing/data/GWAS --ld-pop EUR --outdir results/ --slurm
```
