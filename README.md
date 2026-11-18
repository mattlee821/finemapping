# GWASfinemapping

GWASfinemapping performs locus-level statistical fine-mapping and LD-based variant prioritisation on standardised GRCh38 GWAS-VCFs (specification 1.0.0, as written by [GWASprocessing](https://github.com/mattlee821/GWASprocessing)). It screens each GWAS for genome-wide significant associations ($P \le 5 \times 10^{-8}$), defines merged 3 Mb lead windows, protein `cis_gene` windows and a single extended MHC region, lifts and aligns variants to GRCh37 against the 6 public Pan-UKB LD panels (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`), extracts each locus LD matrix once, and runs up to 5 methods (`p-ld`, `susie-rss`, `susie-rrss`, `carma`, `polyfun`). Multi-gigabyte scratch LD matrices are deleted as each locus finishes; only compact, `coloc`-ready `.rda` (or `.rds`) files, `{metadata,provenance,qc}.json` and a self-contained `report.html` are kept (~25–40 MB per active GWAS). It runs under Nextflow on a Mac or a SLURM cluster, and 2 setup commands install its tools and reference data.

Supporting design notes and literature audits are in [docs/architecture.md](docs/architecture.md), [docs/literature_review.md](docs/literature_review.md), [docs/panukb-ld.md](docs/panukb-ld.md), and [docs/reference-data.md](docs/reference-data.md).

## Installation

You need `bash`, `curl`, and `git` (on a Mac, install the Command Line Tools with `xcode-select --install`), network access to GitHub, Ensembl, UCSC and the public Pan-UKB AWS S3 bucket (`s3://pan-ukb-us-east-1`), and about 20 GB of free disk (~8 GB for the local conda/R runtime, ~3.8 GB for the installed reference data, plus ~5–10 GB of temporary scratch space while loci run).

#### 1. Clone the repository

```bash
git clone https://github.com/mattlee821/finemapping.git GWASfinemapping
cd GWASfinemapping
```

#### 2. Choose local tools or a container

Both routes use the same repository, reference data, manifests and `scripts/finemap.sh` commands; they differ only in where the runtime tools come from. Choose 1 and use its command in step 4.

| Route | Use it when | What step 4 installs |
|---|---|---|
| Local tools (default) | on a Mac, or on a Linux cluster without Apptainer | OpenJDK 21, Nextflow, Python 3.12 (`hail`, `ldcov[s3]`, `lazybgen`, `pyarrow`, `pyliftover`), `bcftools`, `gsl`, R 4.6, pinned `susieR` 2.0.0, pinned `CARMA`, and the local `gwasfinemapping` R/C package under `~/.cache/gwasfinemapping/` (or `GWASFM_TOOLS_ROOT`); no administrator rights needed |
| Container | on a cluster with Apptainer (`module load apptainer` first if your cluster uses modules), or on a machine with Docker | the published image `ghcr.io/mattlee821/finemapping:latest`, inside which tasks execute, plus the host runtime used to launch Nextflow |

#### 3. Choose where the data go (optional)

Tools, reference lookups, Nextflow work directories and caches live outside the repository by default (`~/.cache/gwasfinemapping`), so keeping the repository inside a cloud-synced folder will not sync multi-gigabyte work files. To place them on a dedicated scratch or project volume on a server, copy `.env.example` to `.env` before step 4 (or edit `.env` afterwards):

| Setting in `.env` | Holds |
|---|---|
| `GWASFM_WORK_ROOT` | the fine-mapping outputs, per-GWAS `logs/`, and `manifest.html` (default `<repo>/finemapping`, or `~/GWASfinemapping-data/finemapping` when the repository is inside a cloud-synced folder) |
| `GWASFM_TOOLS_ROOT` | the private micromamba environment (`env/`), `renv-library/`, and pulled `.sif` image (default `~/.cache/gwasfinemapping`) |
| `GWASFM_REFERENCE_ROOT` | the installed reference data (~3.8 GB: chain, GRCh37 FASTA, GRCh38 gene map, 6 Pan-UKB whole-population variant indexes, and PolyFun UKB functional priors; default `<GWASFM_TOOLS_ROOT>/reference`) |
| `GWASFM_CACHE_DIR` | per-run logs (`logs/runs/<run_id>/`), regional index cuts (`panukb/`), and temporary per-locus scratch LD matrices (`scratch/<session_id>/`, deleted automatically as each locus finishes) |
| `GWASFM_NXF_WORK` | Nextflow task work directories and LevelDB cache (default `<GWASFM_TOOLS_ROOT>/nextflow-work`) |
| `GWASFM_SLURM_PARTITION`, `GWASFM_SLURM_ACCOUNT` | optional SLURM `--partition` and `--account` passed to worker jobs under `--slurm` / `-profile slurm` |
| `GWASFM_MAX_FORKS` | maximum concurrent tasks per process / SLURM queue size (default `4` in `.env.example`; override per run with `--max-forks`) |

Keep `GWASFM_REFERENCE_ROOT`, `GWASFM_CACHE_DIR` and `GWASFM_NXF_WORK` on fast local or cluster scratch storage that compute nodes can reach, never on cloud-synced folders. `.env` is ignored by git.

#### 4. Run setup

Install the runtime environment and R packages:

```bash
bash scripts/source.sh                        # local tools + R packages (susieR 2.0.0, CARMA, gwasfinemapping)
bash scripts/source.sh --container apptainer  # container route, with Apptainer
bash scripts/source.sh --container docker     # container route, with Docker
```

Then install and verify the reference data (~3.8 GB total):

```bash
# Single machine or login node (builds all 6 Pan-UKB indexes sequentially; ~2.5 hours first time):
bash scripts/install_reference_data.sh install
bash scripts/install_reference_data.sh check

# Faster on a SLURM cluster (~26 minutes total by building the 6 Pan-UKB indexes in parallel):
bash scripts/install_reference_data.sh --only genome,annotation,polyfun
for p in AFR AMR CSA EAS EUR MID; do
  sbatch --job-name="gwasfm_${p}" scripts/install_reference_data.sbatch --only panukb --populations "${p}"
done
# Once the 6 jobs finish, verify everything is in place:
bash scripts/install_reference_data.sh check
```

Both setup scripts are idempotent: rerunning them skips anything whose checksum or stamp already matches. Individual parts of `scripts/source.sh` can also be run alone:

```bash
bash scripts/source.sh --env-only        # only create .env from .env.example if absent
bash scripts/source.sh --tools-only      # only install or update the conda/Python runtime
bash scripts/source.sh --r-packages      # only restore renv.lock, CARMA, and gwasfinemapping
bash scripts/source.sh --reference-only  # only install/verify reference data
bash scripts/source.sh --force           # rebuild the conda environment and re-pull container
source scripts/source.sh                 # activate the installed tools in your current shell
```

#### 5. Run the tests

About 3–4 minutes; every stage must pass (`=== ALL TESTS PASSED ===`). The suite runs:
1. `[1/4]` 45 Python unit tests (liftover, allele-swap sign invariance, Kanai 3 Mb window merging, single-region MHC policy, PyArrow regional index cuts, manifest preflight and reference checks).
2. `[2/4]` `testthat` unit tests across all R modules (`input`, `prioritisation`, `susie`, `carma-polyfun`).
3. `[3/4]` `R CMD check --no-manual` on the compiled `gwasfinemapping` package (`0 errors | 0 warnings | 0 notes`).
4. `[4/4]` End-to-end offline Nextflow runs testing batch `--manifest` execution across all 5 methods, idempotent state skipping, `--exclude-mhc`, deferred `--report false` + `render_report.R`, user-selected method subsets (`--methods "p-ld susie-rrss polyfun" --output-format rds`), per-locus scratch cleanup, and `coloc.susie` matrix naming.

```bash
bash scripts/test.sh
```

#### 6. Run the benchmark

`scripts/benchmark.sh` validates the installation end-to-end against stored expected summaries in `benchmark/expected/`:

```bash
# Offline zero-network synthetic benchmark (~1 minute; runs all 5 methods):
bash scripts/benchmark.sh --offline

# Real-data live Pan-UKB S3 benchmark (Locke 2015 BMI + Sun 2023 UKB-PPP NPPB P16860_OID20049, --max-loci 2):
bash scripts/benchmark.sh
```

For the real-data benchmark, `scripts/benchmark.sh` automatically looks for standardised `locke_2015_25673413/standardise/GWAScatalog_GCST002783` and `sun_2023_37794186/standardise/P16860_OID20049` under `~/GWASprocessing/GWAS`, `~/research-local/GWASprocessing/data/GWAS`, or `/Volumes/drive1/research/GWASprocessing/data/GWAS` (or pass `--gwas-root /path/to/GWASprocessing/GWAS`), downloads the Pan-UKB EUR LD blocks live from S3, runs all 5 methods, and checks variant-by-variant PIP and credible-set concordance (`Result: PASS`). On a server you can also run the 3-point PyArrow-cut vs Hail parity check on any standardised GWAS-VCF:

```bash
bash scripts/benchmark.sh --server-panukb /home/leem/GWASprocessing/GWAS/sun_2023_37794186/standardise/P16860_OID20049/P16860_OID20049.gwas.vcf.gz
```

## Usage

#### 1. Supply GWAS inputs

`scripts/finemap.sh` reads `GWASprocessing` output folders directly — no file renaming or intermediate conversion is needed. Supply `--outdir` plus 1 of `--input-dir`, `--manifest`, or `--vcf`:

1. **Directory discovery (`--input-dir`):** Point at any of three levels in a `GWASprocessing` tree:
   - the `GWASprocessing` work root (`/home/leem/GWASprocessing/GWAS`) to discover every `<STUDYID>/standardise/<GWAS_ID>/` across all studies (optionally filtered with `--only-study <STUDYID>`);
   - a single study directory (`/home/leem/GWASprocessing/GWAS/sun_2023_37794186` or `.../sun_2023_37794186/standardise`) to discover all GWAS in that study (e.g. all UKB-PPP proteins);
   - a single GWAS folder (`/home/leem/GWASprocessing/GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783`) to run 1 trait.
   Both `standardised.gwas.vcf.gz` and `<output_name>.gwas.vcf.gz` (plus `.tbi`, `.metadata.json` and `.qc.json`) are recognised automatically.
2. **Manifest (`--manifest`):** A tab-separated file with 1 row per GWAS, useful when combining traits from different locations or overriding per-row settings:

   ```text
   study_id	gwas_id	gwas_dir	ld_population	gene_id
   locke_2015_25673413	GWAScatalog_GCST002783	/home/leem/GWASprocessing/GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783	EUR	
   sun_2023_37794186	P16860_OID20049	/home/leem/GWASprocessing/GWAS/sun_2023_37794186/standardise/P16860_OID20049	EUR	P16860
   ```

   | Column | Required | Meaning |
   |---|---|---|
   | `gwas_dir` or `vcf` | yes (one of the two) | path to a `GWASprocessing` `standardise/<GWAS_ID>` folder (`gwas_dir`) or directly to a `.gwas.vcf.gz` (`vcf`) |
   | `study_id`, `gwas_id` | optional with `gwas_dir`; recommended with `vcf` | output hierarchy `<outdir>/<study_id>/<gwas_id>/` (inferred from `metadata.json` or parent folders when omitted) |
   | `ld_population` | no | per-row Pan-UKB panel override (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`; defaults to `--ld-pop`) |
   | `gene_id` | no | UniProt accession (`P16860`) or HGNC symbol (`NPPB`) for `cis_gene` windowing (auto-resolved from `gwas_id` or `metadata.json` for `UNIPROT_OID` / `UNIPROT_SEQID` protein IDs) |
   | `metadata`, `qc`, `trait`, `ancestry` | no | explicit sidecar paths or metadata overrides when not alongside `vcf` |

3. **Single VCF (`--vcf`):** Pass `--vcf /path/to/trait.gwas.vcf.gz` with optional `--study-id`, `--gwas-id`, `--metadata`, `--qc`, `--ancestry`, and `--gene-id`.

#### 2. Choose fine-mapping methods and options

All 5 methods run by default (`--methods "p-ld susie-rss susie-rrss carma polyfun"`). Pass a space- or comma-separated list to run only a subset.

| Option | Default | Meaning |
|---|---|---|
| `--methods LIST` | `"p-ld susie-rss susie-rrss carma polyfun"` | methods to run (`p-ld`, `susie-rss`, `susie-rrss`, `carma`, `polyfun`; hyphens or underscores accepted) |
| `--ld-pop POP` | `EUR` | Pan-UKB LD reference panel: `AFR`, `AMR`, `CSA`, `EAS`, `EUR`, or `MID` |
| `--allow-ancestry-mismatch` | `false` | allow a single-ancestry GWAS (e.g. `EAS`) to run against a different single-ancestry LD panel (multi-ancestry `ALL` GWAS are allowed automatically and flagged in `qc.json` and `report.html`) |
| `--exclude-mhc` | `false` | skip the extended MHC (`chr6:25,500,000–34,000,000` GRCh38) instead of fine-mapping it as 1 merged region |
| `--significance-lp FLOAT` | `7.30103` ($P \le 5 \times 10^{-8}$) | $-\log_{10} P$ threshold required to activate any region (`lead_locus`, `cis_gene`, or `mhc`) |
| `--cis-significance-lp FLOAT` | unset (inherits `--significance-lp`) | optional $-\log_{10} P$ override for `cis_gene` windows only, recorded in `provenance.json` |
| `--region-flank-bp INT` | `1500000` | half-width around greedy lead variants in bp (`±1.5 Mb` = 3 Mb nominal window, merged on overlap up to 9.5 Mb) |
| `--max-loci INT` | `0` (all) | cap the number of candidate regions analyzed per GWAS (`cis_gene` ordered first, then descending lead `LP`) |
| `--output-format FMT` | `rda` | saved method container format: `rda`, `rds`, or `both` |
| `--compress COMP` | `xz` | compression for `.rda`/`.rds` method files: `xz` (~35% smaller at 10,000-GWAS scale) or `gzip` (faster write) |
| `--row-limit INT` | `0` (all) | run at most `INT` incomplete GWAS rows from the discovered batch |
| `--max-forks INT` | `$GWASFM_MAX_FORKS` or `2` | maximum concurrent tasks / SLURM queue size |
| `--report true\|false` | `true` | render self-contained per-GWAS `report.html` during finalisation |
| `--keep-intermediates` | `false` | retain per-locus lifted VCFs, mapping tables, and `matrix.npy` files under `<outdir>/<study_id>/<gwas_id>/loci/` |
| `--dry-run` | `false` | run manifest discovery and preflight validation without launching fine-mapping |
| `--force` | `false` | re-run rows even when complete outputs exist for the current settings hash |

#### 3. Run it

When `--outdir` is omitted, outputs default to `$GWASFM_WORK_ROOT` (`~/GWASfinemapping/finemapping/<STUDYID>/<GWASID>/`).

```bash
# Dry-run preflight on a single study from GWASprocessing
bash scripts/finemap.sh \
  --input-dir /home/leem/GWASprocessing/GWAS/sun_2023_37794186 \
  --ld-pop EUR \
  --dry-run

# Local execution: 1 GWAS, all 5 methods
bash scripts/finemap.sh \
  --input-dir /home/leem/GWASprocessing/GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783 \
  --ld-pop EUR

# Local or interactive node: custom subset of methods
bash scripts/finemap.sh \
  --input-dir /home/leem/GWASprocessing/GWAS/sun_2023_37794186 \
  --methods "p-ld susie-rrss polyfun" \
  --ld-pop EUR

# Submit orchestrator to SLURM (fans out SCREEN_GWAS, PREP_LOCUS, and method tasks as SLURM jobs)
sbatch scripts/finemap.sh \
  --manifest manifests/sun_mr18.tsv \
  --ld-pop EUR \
  --max-forks 16 \
  --slurm
```

#### 4. Check what has been processed

Open `<outdir>/manifest.html` (by default `~/GWASfinemapping/finemapping/manifest.html`). Updated after every run, it summarises every GWAS in `<outdir>` with its study, trait, LD panel, status badge, analyzed region count, 95% credible-set counts, and links to per-GWAS reports. Each GWAS writes its outputs and logs to `<outdir>/<STUDYID>/<GWAS_ID>/`:

```text
finemapping/
  manifest.html                       batch dashboard across all studies in <outdir>
  .finemap_state.tsv                  idempotent completion ledger (keyed by settings_hash)
  <STUDYID>/
    <GWAS_ID>/
      metadata.json                   upstream study metadata + fine-mapping summary counts
      provenance.json                 build, chain SHA-256, Pan-UKB panel N, thresholds, runtime
      qc.json                         per-locus mapping/LD/method diagnostics + step_telemetry
      report.html                     self-contained HTML report (tables, flags, base64 locus plots)
      logs/                           01_screen_gwas.log, 02_prep_*.log, 03_method_*.log,
                                      04_finalise_gwas.log, trace.tsv, nextflow.*.log, slurm-*.out
      methods/
        p-ld.rda                      greedy r²-pruned sets at r² < 0.1, 0.01, 0.001
        susie-rss.rda                 SuSiE-RSS PIPs, 95% CSs, diagnostics, coloc-ready fit
        susie-rrss.rda                robust SuSiE-rRSS PIPs, 95% CSs, eb_mix mismatch diagnostics
        carma.rda                     CARMA shotgun PIPs, 95% CSs, Bayes-factor outlier flags
        polyfun.rda                   PolyFun + SuSiE functionally weighted PIPs, 95% CSs, fit
```

Only the method files requested via `--methods` (in `.rda`, `.rds`, or both per `--output-format`) are written under `methods/`.

#### 5. Rerun when something changes

`scripts/prepare_manifest.py` hashes each GWAS-VCF size, `ld_population`, `gene_id`, requested `methods`, `output_format`, significance thresholds, window parameters, MHC policy, `L`, and `max_iter` into a 16-character `settings_hash` stored in `<outdir>/.finemap_state.tsv`. Rerunning the same command skips already-complete GWAS in under 1 second without starting the Nextflow JVM. Changing `--methods` (for example adding `carma` or `polyfun`), changing `--ld-pop`, or passing `--force` automatically re-queues the affected GWAS. A failure in 1 GWAS is isolated to its `.row_record.json` and reported in the end-of-run summary without aborting sibling GWAS in the batch.

If a large batch was run with `--report false`, render any individual GWAS report later purely from its `.rda`/`.rds` and `.json` files:

```bash
source scripts/source.sh
Rscript scripts/render_report.R --gwas-dir /home/leem/GWASfinemapping/results/sun_2023_37794186/P16860_OID20049
```

#### 6. Read results in R, run `coloc.susie`, or step through scripts interactively

Every Bayesian method file (`susie-rss.rda`, `susie-rrss.rda`, `polyfun.rda`) names `colnames(fit$alpha)`, `colnames(fit$lbf_variable)` and `names(fit$pip)` with canonical GRCh38 `chr:pos:ref:alt` (`GWAS_VARIANT_ID`), so `coloc::coloc.susie()` and `coloc::coloc.bf_bf()` work directly across traits without needing an LD matrix:

```r
e_pqtl <- new.env(parent = emptyenv())
e_crc  <- new.env(parent = emptyenv())
load("results/sun_2023_37794186/P16860_OID20049/methods/polyfun.rda", envir = e_pqtl)
load("results/huyghe_2021_33632709/colon_overall/methods/polyfun.rda", envir = e_crc)

fit_pqtl <- e_pqtl$polyfun$regions[["cis_NPPB_chr1_10857464_12858945"]]$result$fit
fit_crc  <- e_crc$polyfun$regions[[1]]$result$fit

coloc_res <- coloc::coloc.susie(fit_pqtl, fit_crc)
```

To inspect or debug locus preparation and each fine-mapping method line-by-line in RStudio or Positron, open:
- `scripts/interactive_prep_locus.R` (steps through Stage 1 screening, Stage 2 liftover/LD extraction, and inspects `matched_variants.tsv` and `matrix.npy`)
- `scripts/methods/run_p_ld.R`
- `scripts/methods/run_susie_rss.R`
- `scripts/methods/run_susie_rrss.R`
- `scripts/methods/run_carma.R`
- `scripts/methods/run_polyfun.R`

## Pipeline steps

| Stage | Script / Process | What happens | Reference / Method |
|---|---|---|---|
| 1. Screen | `scripts/screen_gwas.py` (`SCREEN_GWAS`) | Stream variants at $P \le 5 \times 10^{-8}$; exit early if 0 hits; build greedy 3 Mb merged lead windows, `cis_gene` (`gene ± 1 Mb`), and single-region MHC (`chr6:25.5–34.0 Mb`) | Kanai et al. 2021; Sun et al. 2023 |
| 2. Prepare locus | `scripts/prep_locus.py` (`PREP_LOCUS`) | Query region from GWAS-VCF, lift GRCh38→GRCh37, left-align/verify alleles against GRCh37 FASTA, cut regional Pan-UKB index with PyArrow, extract & orient Pan-UKB BlockMatrix LD into Fortran-order `matrix.npy` | `pyliftover`, `bcftools norm`, `ldcov` |
| 3a. `p-ld` | `scripts/methods/run_p_ld.R` (`RUN_P_LD`) | Deterministic greedy $r^2$ pruning at $r^2 < 0.1, 0.01, 0.001$ from shared `matrix.npy` | Selected-Push Sieve (`R/prioritize_pld.R`) |
| 3b. `susie-rss` | `scripts/methods/run_susie_rss.R` (`RUN_BAYESIAN_METHOD`) | Summary-statistic Sum of Single Effects ($L = 10$, `max_iter = 50`, 95% CS with `min_abs_corr >= 0.5`) using in-place C $X^\top X$ scaling | Zou et al. 2022 (`susieR` 2.0.0) |
| 3c. `susie-rrss` | `scripts/methods/run_susie_rrss.R` (`RUN_BAYESIAN_METHOD`) | Finite-reference (`R_finite`) and empirical-Bayes mixture LD-mismatch (`R_mismatch = "eb_mix"`) SuSiE-rRSS | `susieR` 2.0.0 |
| 3d. `carma` | `scripts/methods/run_carma.R` (`RUN_BAYESIAN_METHOD`) | Shotgun stochastic search spike-and-slab fine-mapping with Bayes-factor LD-discrepancy outlier detection (`outlier.BF.index = 1/3.2`) | Yang et al. 2023 (`ZikunY/CARMA` + C/LAPACK batch evaluator) |
| 3e. `polyfun` | `scripts/methods/run_polyfun.R` (`RUN_BAYESIAN_METHOD`) | Functionally informed SuSiE-RSS weighted by UK Biobank Baseline-LF2.2 meta-analyzed per-SNP heritability priors (`snpvar_meta`, $q = 100$ floor) | Weissbrod et al. 2020 (`omerwe/polyfun` + `susieR` 2.0.0) |
| Clean scratch | `CLEAN_LOCUS_SCRATCH` | Unlink `matrix.npy` and `matrix.ids.txt` immediately after all requested methods finish for a locus (unless `--keep-intermediates`) | Nextflow head / local cleanup |
| 4. Finalise | `scripts/finalise_gwas.py` (`FINALISE_GWAS`) | Combine per-locus `.rds` results into GWAS-level `.rda`/`.rds` files in parallel, reconcile telemetry into `{metadata,provenance,qc}.json`, render `report.html`, update `.finemap_state.tsv` and `manifest.html` | `aggregate_method_results.R`, `render_report.R` |

### 1. Screen and window definition

`scripts/screen_gwas.py` scans the indexed GRCh38 GWAS-VCF with 4-thread `bcftools view` for variants meeting `--significance-lp` (default `7.30103`, $P \le 5 \times 10^{-8}$). When a GWAS has 0 significant variants (and no override on `cis_gene`), the pipeline writes valid empty method containers and `overall_status = "NO_SIGNIFICANT_LOCI"` in ~30 seconds without touching liftover, S3, or R.

Active associations are partitioned into three region types:
- **`cis_gene` (`gene_start - 1 Mb` to `gene_end + 1 Mb`):** For protein GWAS (`UNIPROT_OLINKID` or `UNIPROT_SEQID`, with UniProt isoform suffixes such as `O43521-2` stripped before lookup, or explicit `--gene-id`), coordinates are looked up in `gene_regions_GRCh38.tsv` (11,105 unambiguous human UniProt→GRCh38 mappings; multi-gene protein complexes listed in `gene_regions_GRCh38.ambiguous.tsv` skip `cis_gene` and fall back cleanly to `lead_locus`). Activated only when $\max(\text{LP})$ in the cis window meets `--cis-significance-lp` (default $P \le 5 \times 10^{-8}$).
- **`mhc` (`chr6:25,500,000–34,000,000` GRCh38):** Following Sun et al. (2023), any significant association in the extended MHC activates a single 8.50 Mb locus (~8.47 Mb on GRCh37, safely inside Pan-UKB's 10 Mb stored LD band) rather than fragmenting long-range HLA LD across overlapping windows. Every variant carries `mhc_region` (extended `25.5–34.0 Mb`) and `mhc_core_grc` (NCBI GRC core MHC `chr6:28,510,120–33,480,577` GRCh38 / `chr6:28,477,797–33,448,354` GRCh37) boolean flags. Pass `--exclude-mhc` to skip it.
- **`lead_locus` (Kanai greedy `±1.5 Mb` = 3 Mb merged):** Remaining genome-wide significant variants outside `mhc` are ranked by descending `LP`, placed in `±1.5 Mb` windows, and merged when overlapping. Any merged chain wider than `--max-region-width-bp` (`9,500,000` bp) is split at its lowest-`LP` valley so all pairwise correlations remain inside Pan-UKB's 10 Mb BlockMatrix radius.

### 2. Shared locus preparation and Pan-UKB LD extraction

For each selected locus, `scripts/prep_locus.py` queries all GWAS-VCF records in the GRCh38 span (preserving `ES`, `SE`, `LP`, `SS`, `AF`, `AFS`, and `INFO/HM`), lifts every nucleotide of `REF` to GRCh37 via `pyliftover` (complementing reverse-strand chain blocks), and runs `bcftools norm --check-ref s` against the Ensembl GRCh37 primary assembly FASTA so indels are left-aligned and reference mismatches are caught and counted rather than mis-mapped.

Lifted GRCh37 variants are matched against a regional PyArrow slice cut from the installed whole-population Pan-UKB variant index (`UKBB.<POP>.variant_index.parquet`, taking milliseconds offline without booting Spark/Hail). Both exact `(REF, ALT)` and allele-swapped `(ALT, REF)` panel matches are resolved: when a study's GRCh38 `REF/ALT` is reversed relative to the Pan-UKB GRCh37 entry (`PANEL_ALLELE_SWAPPED = TRUE`), the extracted LD row and column signs are flipped while `ES` and `EZ` stay oriented to the GWAS `ALT` allele. Duplicate GWAS variants mapping to the same `PANEL_IDX` keep the highest-`LP` record and log excluded counts in `qc.json`. The submatrix is fetched from `s3://pan-ukb-us-east-1/ld_release/UKBB.<POP>.ldadj.bm` (8 parallel workers), symmetrised and unit-diagonal normalised in place, and written directly as a Fortran-order binary float64 `matrix.npy` shared by all downstream method scripts.

### 3. Prioritisation and fine-mapping methods

Each method runs in its own isolated R process (`R_PROFILE_USER=""` with pre-resolved `R_LIBS_USER` for <0.1s startup) reading the shared `matched_variants.tsv` and `matrix.npy`:
- **`p-ld` (`run_p_ld.R`):** Ranks variants by descending `LP` (tie-broken deterministically by coordinate and `GWAS_VARIANT_ID`) and greedily selects variants with $r^2 < t$ against all previously selected variants independently at $t \in \{0.1, 0.01, 0.001\}$. Reports heuristic prioritisation sets only (never PIPs or Bayesian credible sets).
- **`susie-rss` (`run_susie_rss.R`):** Scales `R` to $X^\top X = (n_{\text{eff}} - 1)R$ in place in C (`src/fast_matrix.c`, avoiding 4 full-matrix copies inside `susieR`), fits `susieR` 2.0.0 (`L = 10`, `max_iter = 50`, $n_{\text{eff}} = \min(\text{SS})$), and filters 95% credible sets at exact submatrix purity `min_abs_corr >= 0.5`.
- **`susie-rrss` (`run_susie_rrss.R`):** Fits mismatch-aware `SuSiE-rRSS` in `susieR` 2.0.0 with finite-reference adjustment (`R_finite` set to the exact Pan-UKB panel sample size) and empirical-Bayes residual mixture modeling (`R_mismatch = "eb_mix"`), recording `SER_fallback` / `R_reliability_flag` and `artifact_flag` diagnostics.
- **`carma` (`run_carma.R`):** Runs shotgun stochastic search Bayesian fine-mapping (`ZikunY/CARMA`, default `Spike-slab` effect prior, Poisson model-size prior, `L = 10`, deterministic RNG seed, `outlier.switch = TRUE`, `outlier.BF.index = 1/3.2`) accelerated by a vectorised LAPACK `dsyev` batch marginal-likelihood evaluator in `src/fast_matrix.c` and closed-form sparse priors, recording flagged outlier variant IDs and indices.
- **`polyfun` (`run_polyfun.R`):** Matches lifted GRCh37 `(CHR, BP, REF, ALT)` variants (either allele order, autosomes 1–22) against the 19.48 million precomputed UK Biobank Baseline-LF2.2 meta-analyzed per-SNP heritabilities (`snpvar_meta`, Weissbrod et al. 2020 Approach 1, $q = 100$ floor $\max_j \sigma_j^2 / 100$; unmatched rare variants receive the locus baseline floor and are counted in `qc.json`) before loading `matrix.npy`, then fits prior-weighted `SuSiE-RSS` (`prior_weights = snpvar / sum(snpvar)`).

### 4. Aggregation, provenance, QC and reporting

`scripts/finalise_gwas.py` aggregates all per-locus `region_result.rds` files into GWAS-level `.rda` and/or `.rds` method files in parallel, collects per-step wall/CPU/memory metrics into `qc.json` (`step_telemetry`), assigns `overall_status` (`PASS`, `WARN_DIAGNOSTIC_FLAGS` when any locus fails to converge or triggers `susie-rrss` / `carma` LD-mismatch flags, `NO_SIGNIFICANT_LOCI`, or `NO_PANEL_MATCHED_LOCI`), and renders `report.html` (summary cards, region/method/flag table, 95% credible-set lead table across all active Bayesian methods, and 3-panel base64 PNG plots per locus).

## Output schemas

Inside each Bayesian `.rda` (`load("susie-rss.rda")` loads `susie_rss`; `susie-rrss.rda` → `susie_rrss`; `carma.rda` → `carma`; `polyfun.rda` → `polyfun`) or `.rds` (`readRDS(...)`):

- `<obj>$method_id`: `"susie_rss"`, `"susie_rrss"`, `"carma"`, or `"polyfun"`.
- `<obj>$regions[[region_id]]`: list with `region_id`, `region_type` (`"cis_gene"`, `"lead_locus"`, or `"mhc"`), and `result`:
  - `result$variants`: data frame containing GRCh38 `GWAS_VARIANT_ID`, `GWAS_CHROM`, `GWAS_POS`, `GWAS_REF`, `GWAS_ALT`, `ID` (rsID), `LP`, `ES`, `SE`, `SS`, `AF`, `AFS`, `HM`, `EZ`, lifted GRCh37 `LD_CHROM`, `LD_POS`, `LD_REF`, `LD_ALT`, `PANEL_IDX`, `PANEL_ALLELE_SWAPPED`, `LD_R2_TO_LEAD`, `cis_trans`, `mhc_region`, `mhc_core_grc`, `PIP`, `credible_set_id` ( Plus `OUTLIER` in `carma`; `SNPVAR`, `PRIOR_WEIGHT`, `PRIOR_MATCHED` in `polyfun`).
  - `result$credible_sets`: 1 row per purity-passing 95% credible set (`credible_set_id`, `effect_index`, `size`, `coverage`, `min_abs_corr`, `mean_abs_corr`, `median_abs_corr`).
  - `result$credible_set_membership`: long-format mapping of `credible_set_id` to `GWAS_VARIANT_ID`.
  - `result$diagnostics`: convergence, iteration count, ELBO / Bayes-factor summaries, sample-size heterogeneity warnings, and mismatch/outlier/prior-coverage metrics.
  - `result$fit` (`susie_rss`, `susie_rrss`, `polyfun`): `susie` S3 object (`O(L * M)`, LD matrix omitted) with `GWAS_VARIANT_ID` column/element names on `alpha`, `lbf_variable`, `mu`, `mu2`, and `pip` for direct use with `coloc::coloc.susie()`.

Inside `p-ld.rda` (`p_ld$regions[[region_id]]$result`): long-format data frame stacked across `r2_threshold` $\in \{0.1, 0.01, 0.001\}$ with `rank`, `selected`, `pruned_by`, `max_r2_to_selected`, and all GWAS/liftover variant columns.

## Reference data

Installed once under `GWASFM_REFERENCE_ROOT` by `bash scripts/install_reference_data.sh` and checksum-verified by `reference_manifest.json`:

| Under `GWASFM_REFERENCE_ROOT` | Size | Use | Source |
|---|---|---|---|
| `genome/hg38ToHg19.over.chain.gz` | 1.2 MB | GRCh38→GRCh37 coordinate liftover | UCSC |
| `genome/GRCh37.fa`, `GRCh37.fa.fai` | 860 MB | GRCh37 reference verification and indel left-alignment | Ensembl release 113 primary assembly |
| `annotation/gene_regions_GRCh38.tsv` (and `.ambiguous.tsv`) | 0.9 MB | UniProt / HGNC `cis_gene` coordinates on GRCh38 (11,105 unambiguous genes) | `usefuldata` `mapping_GRCh38_p14_uniprot` (`30f2e833`) |
| `panukb/UKBB.<POP>.variant_index.parquet` (+ `.json`) | 2.8 GB (all 6 pops) | Offline regional PyArrow cuts for `AFR` ($N = 6{,}637$), `AMR` ($982$), `CSA` ($8{,}876$), `EAS` ($2{,}709$), `EUR` ($420{,}542$), `MID` ($1{,}599$) | Pan-UKB `s3://pan-ukb-us-east-1/ld_release/UKBB.<POP>.ldadj.variant.ht` |
| `polyfun/by_chrom/chr1..22.parquet` (+ `polyfun_meta.json`) | 60.8 MB | 19,476,620 autosomal GRCh37 UKB meta-analyzed functional priors (`snpvar_meta`) | Weissbrod et al. 2020 (`omerwe/polyfun@8b17b4ec`, Alkes Group broadinstitute mirror) |
| *(live on demand)* `s3://pan-ukb-us-east-1/ld_release/UKBB.<POP>.ldadj.bm` | streamed per locus | Covariate-adjusted signed in-sample LD BlockMatrix (10 Mb radius) | Pan-UKB public AWS S3 (`us-east-1`, anonymous access) |
