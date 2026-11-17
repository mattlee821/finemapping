# Pan-UKB LD panels and reference workflow

The Pan-UKB LD BlockMatrices and companion variant-index Hail Tables (`s3://pan-ukb-us-east-1/ld_release/`) are publicly readable from AWS S3 without requester-pays charges.

## Released populations, variant counts, and exact sample sizes

Pan-UKB releases six ancestry-specific GRCh37 LD panels (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`), stored as covariate-adjusted upper-triangular BlockMatrices with a **10 Mb stored radius** (`10,000,000 bp`). There is no `ALL` BlockMatrix: pooled multi-ancestry meta-analyses must not be treated as ancestry-matched to a single Pan-UKB panel without an explicit justification.

Sample sizes (`n_samples`, verified directly from `UKBB.<POP>.ldadj.variant.ht/globals`) and installed biallelic variant counts:

| Population | Pan-UKB `n_samples` (used for `SuSiE-rRSS` `R_finite`) | Biallelic indexed variants (GRCh37) |
|---|---:|---:|
| `AFR` | `6,637` | `22,082,831` |
| `AMR` | `982` | `11,713,690` |
| `CSA` | `8,876` | `16,161,130` |
| `EAS` | `2,709` | `10,225,280` |
| `EUR` | `420,542` | `23,960,350` |
| `MID` | `1,599` | `14,127,512` |

## One-time installation vs per-locus PyArrow slicing

Building regional indexes through Hail + Spark against S3 on every locus is slow and redundant at >10,000-GWAS scale. Instead:

1. **One-time reference install (`scripts/install_reference_data.sh`):** Builds all six whole-population Parquet variant indexes (`panukb/UKBB.<POP>.variant_index.parquet` + `.json`, ~3.6 GB total including GRCh37 FASTA, UCSC `hg38ToHg19.over.chain.gz`, and UniProt GRCh38 gene coordinates) once per machine/cluster under `GWASFM_REFERENCE_ROOT` (default `~/.cache/gwasfinemapping/reference`). See [docs/reference-data.md](reference-data.md).
2. **Fast offline regional cut (`cut_regional_index` in `scripts/run_finemapping.py`):** At runtime, PyArrow slices only the GRCh37 contig/position row groups overlapping a locus in seconds with zero Hail/Spark overhead, writing a small regional Parquet slice for `ldcov` (`s3://pan-ukb-us-east-1/ld_release/UKBB.<POP>.ldadj.bm`). If no whole-population index is installed, the runner falls back to `scripts/build_panukb_variant_index.sh`.

## 10 Mb BlockMatrix radius and MHC / merged-window constraints

Because Pan-UKB BlockMatrices store pairwise correlations only up to **10 Mb** on GRCh37:
- **MHC single region:** Extended HLA is defined on GRCh38 as **`chr6:25,500,000–34,000,000`** (`8.50 Mb`, lifting to ~`chr6:25,499,000–33,968,000` on GRCh37 = `8.47 Mb < 10.0 Mb`), matching Sun et al. (2023) while keeping every variant pair strictly inside the 10 Mb stored band. Adjacent non-MHC windows are truncated at `25,499,999` / `34,000,001` so chain-merging cannot push MHC above 10 Mb.
- **Hyper-polygenic merged loci:** Non-MHC greedy `±1.5 Mb` (3 Mb) merged windows are guarded by `--max-region-width-bp 9500000` (splitting rare >9.5 Mb chains at the deepest $-\log_{10} P$ valley between leads).
- **Missing LD policy:** Out-of-band LD is unstored, never zero; `--on-missing error` enforces complete within-region LD coverage.

## Storage-light output layout

Each fine-mapped GWAS writes only self-contained `xz`-compressed `methods/{p-value,p-ld,susie-rss,susie-rrss}.rda`, `{metadata,provenance,qc}.json`, and `report.html` (~20–24 MB per GWAS), deleting temporary `.npz`/`.npy` LD matrices and lifted VCFs immediately after each locus completes (unless `--keep-intermediates` is set). See [docs/architecture.md](architecture.md) for the complete output and `coloc.susie` contract.
