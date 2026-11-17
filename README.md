# GWASfinemapping

Scalable statistical fine-mapping and variant prioritisation for standardised GRCh38 GWAS-VCF summary statistics (designed for both individual complex-trait GWAS and >10,000-trait proteogenomic studies) against the six public Pan-UKB LD panels (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`).

`GWASfinemapping` consumes the standardised GWAS-VCF outputs of `GWASprocessing` without runtime coupling: see [docs/architecture.md](docs/architecture.md), [docs/literature_review.md](docs/literature_review.md), [docs/panukb-ld.md](docs/panukb-ld.md), and [docs/reference-data.md](docs/reference-data.md).

## Key features

- **Evidence-backed region & MHC defaults:** Default $P \le 5 \times 10^{-8}$ (`LP >= 7.30103`) significance gate across all region types; Kanai et al. (2021) greedy `±1.5 Mb` (**3 Mb**) merged lead windows (with a 9.5 Mb ceiling respecting Pan-UKB's 10 Mb BlockMatrix radius); gene-anchored `cis_gene` (`gene ± 1 Mb`) windows for molecular traits; and single-region extended MHC (`chr6:25.5–34.0 Mb` GRCh38) with `mhc_region` / `mhc_core_grc` flags and `--exclude-mhc`.
- **Ordered Phase 1 methods:** Deterministic `p-value` ranking and `p-ld` ($r^2 < 0.1, 0.01, 0.001$) prioritisation heuristics, followed by Bayesian `SuSiE-RSS` and mismatch-aware `SuSiE-rRSS` (`susieR` 2.0.0, `R_finite`, `R_mismatch = "eb_mix"`).
- **Storage-light, `coloc`-ready outputs (~20 MB/GWAS; ~200 GB per 10,000 GWAS):** Retains only four standalone `xz`-compressed `.rda` files (`methods/{p-value,p-ld,susie-rss,susie-rrss}.rda` with `GWAS_VARIANT_ID`-named `alpha` and `lbf_variable` matrices ready for `coloc::coloc.susie()` and `coloc::coloc.bf_bf()`), `{metadata,provenance,qc}.json`, and a self-contained `report.html` (plus batch `manifest.html`), automatically cleaning up multi-gigabyte scratch LD files.
- **Portable local, container & SLURM setup:** One-command bootstrap (`scripts/source.sh`), reference installer (`scripts/install_reference_data.sh`), validation test suite (`scripts/test.sh`), and deterministic benchmark (`scripts/benchmark.sh`) matching `GWASprocessing`.

## Quickstart

```bash
# 1. Install pinned conda/Python/R environment (and pull Apptainer .sif on HPC)
bash scripts/source.sh
source scripts/source.sh

# 2. Install & verify reference data once (~3.6 GB under ~/.cache/gwasfinemapping/reference)
bash scripts/install_reference_data.sh install
bash scripts/install_reference_data.sh check

# 3. Validate installation (unit tests + end-to-end synthetic smoke test & benchmark)
bash scripts/test.sh
bash scripts/benchmark.sh

# 4. Run fine-mapping locally or on SLURM (single GWAS-VCF or batch manifest)
bash scripts/finemap.sh --vcf /path/to/trait.vcf.gz --study-id study1 --gwas-id trait1 --ld-pop EUR --outdir results/
# Or submit directly to SLURM:
sbatch scripts/finemap.sh --manifest manifest.tsv --ld-pop EUR --outdir results/ --slurm
```
