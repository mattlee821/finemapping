# Architecture and data contracts

**Status:** Cascading 4-stage pipeline and all five methods (`p-ld`, `susie-rss`, `susie-rrss`, `carma`, `polyfun`) finalised and benchmarked (18 Nov 2026); targets >10,000-GWAS storage and CPU efficiency, shared per-locus LD preparation, interactive standalone per-method R scripts, evidence-backed region/MHC defaults, idempotent batch execution over `GWASprocessing` outputs, self-contained HTML reporting, and container/SLURM portability. Reusable epidemiological and statistical fine-mapping defaults live in `~/.claude/skills/methods/references/fine-mapping.md`. See [the Pan-UKB LD guide](panukb-ld.md) and [the reference-data guide](reference-data.md).

## System boundaries & cascading DAG

`GWASprocessing` produces standardised GRCh38 GWAS-VCF files. `GWASfinemapping` consumes those files without runtime coupling and executes a cascading Nextflow DAG (`finemapping.nf`):

```text
PREPARE_MANIFEST -> SCREEN_GWAS --(0 sig hits)--> FINALISE_GWAS -> UPDATE_STATE
                        |
                   (>=1 sig locus)
                        v
                   PREP_LOCUS (parallel per locus: liftover + Pan-UKB match + Fortran .npy LD)
                        |
          +-------------+-------------+-------------+-------------+
          v             v             v             v             v
       run_p_ld   run_susie_rss run_susie_rrss  run_carma    run_polyfun
          +-------------+-------------+-------------+-------------+
                        |
                        v
                   FINALISE_GWAS (combine .rda/.rds, write JSONs + report.html, delete scratch .npy)
                        |
                        v
                   UPDATE_STATE
```

## Components

| Component | Responsibility |
|---|---|
| `scripts/source.sh` + `scripts/lib/env.sh` | Bootstrap micromamba (`Hail`, `ldcov[s3]`, `bcftools`, `gsl`) + compiled R library (`gwasfinemapping`, `susieR` 2.0.0, `CARMA`) + optional Apptainer `.sif`; resolve `.libPaths()` directly to avoid `renv` startup overhead on every method task. |
| `scripts/install_reference_data.sh` | Install and verify UCSC `hg38ToHg19.over.chain.gz`, Ensembl GRCh37 FASTA + `.fai`, UniProt GRCh38 gene coordinates (`usefuldata`), PolyFun UK Biobank per-chromosome functional heritability priors (`polyfun/by_chrom/chr1..22.parquet`), and all 6 Pan-UKB GRCh37 Parquet variant indexes (`AFR`, `AMR`, `CSA`, `EAS`, `EUR`, `MID`). |
| Stage 1: `scripts/screen_gwas.py` | Scan GRCh38 GWAS-VCF at $P \le 5 \times 10^{-8}$ (`LP >= 7.30103` default), build `cis_gene` (`gene ± 1 Mb`), greedy 3 Mb (`±1.5 Mb`) merged genome-wide windows (9.5 Mb ceiling), and single-region extended MHC (`chr6:25.5–34.0 Mb`), exiting early when 0 loci pass. |
| Stage 2: `scripts/prep_locus.py` + `scripts/interactive_prep_locus.R` | Perform regional VCF extraction, GRCh38→GRCh37 liftover + swapped-indel FASTA pre-normalisation + `bcftools norm`, Pan-UKB `VariantIndex` matching, and in-place normalised Fortran-order float64 `matrix.npy` + `matched_variants.tsv` creation shared by all methods. |
| Stage 3: `scripts/methods/run_{p_ld,susie_rss,susie_rrss,carma,polyfun}.R` | Standalone R scripts designed for both interactive line-by-line execution in RStudio/Positron and parallel Nextflow execution (`--methods` accepts any space/comma subset, defaulting to all 5 methods). |
| Stage 4: `scripts/finalise_gwas.py` + `scripts/aggregate_method_results.R` | Combine per-locus `region_result.rds` outputs in parallel into GWAS-level `.rda` and/or `.rds` files (`xz` or `gzip`), record per-step CPU/wall/RSS telemetry in `qc.json`, render `report.html`, and delete scratch `.npy` matrices immediately. |
| R statistical package (`R/` + `src/fast_matrix.c`) | Zero-copy Fortran `.npy` reader with C in-place $X^\top X = (n - 1)R$ scaling (`gwasfm_frob_and_scale_inplace`) and LAPACK batch Spike-slab marginal likelihoods (`gwasfm_carma_spike_slab_batch`); exports `read_gwas_vcf`, `read_ld_matrix`, `validate_gwas_stats`, `validate_ld_matrix`, `prioritize_p`, `prioritize_pld`, `run_susie_rss`, `run_susie_rrss`, `run_carma`, and `run_polyfun_susie`. |

## Storage-light output contract (>10,000 GWAS scale)

Each published `<output_root>/<study_id>/<gwas_id>/` directory retains only:

```text
<study_id>/<gwas_id>/
├── metadata.json          # Study/trait/gene/panel summary & per-method CS counts
├── provenance.json        # Software/git versions, CLI parameters, reference SHA-256s, Pan-UKB sample size
├── qc.json                # Complete locus manifest, liftover/panel attrition, LD validation, CS/outlier flags & step_telemetry
├── report.html            # Self-contained QC & locus visualisation report (base64 plots, zero external assets)
└── methods/
    ├── p-ld.rda           # Standalone xz-compressed R object `p_ld` (r2 = 0.1, 0.01, 0.001)
    ├── susie-rss.rda      # Standalone xz-compressed R object `susie_rss` (coloc.susie-ready `fit`)
    ├── susie-rrss.rda     # Standalone xz-compressed R object `susie_rrss` (R_finite + R_mismatch = "eb_mix")
    ├── carma.rda          # Standalone xz-compressed R object `carma` (PIPs, CSs, CARMA_OUTLIER flags)
    └── polyfun.rda        # Standalone xz-compressed R object `polyfun` (SNPVAR priors + coloc.susie-ready `fit`)
```

*(Only the files corresponding to the requested `--methods` subset are written; `--output-format rds|both` optionally writes `.rds` files.)*

### `coloc.susie` & `coloc.bf_bf` compatibility contract

In `susie-rss.rda`, `susie-rrss.rda`, and `polyfun.rda`, every region's `fit` object sets `colnames(fit$alpha)`, `colnames(fit$lbf_variable)`, `colnames(fit$mu)`, `colnames(fit$mu2)`, and `names(fit$pip)` to canonical GRCh38 `GWAS_VARIANT_ID` (`chr:pos:ref:alt`), alongside `fit$sets` (`min_abs_corr = 0.5`, `coverage = 0.95`), enabling direct `coloc::coloc.susie(fit1, fit2)` execution with no LD matrix required.
