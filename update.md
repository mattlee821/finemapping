# Update for the next agent on finemapping/ (written 17 Nov 2026)

Written for: the agent that works on the finemapping/ project next. It says what was done in the last session, why, and what came out of it, and it hands over one step that has to be redone (the section "Step to redo").

Nothing in this project has been committed. The restructured project (README rewrite, legacy deletions, new scripts, docs and tests) is all uncommitted. The agent working on this project is in charge of committing and pushing it (Matt's decision, 17 Nov 2026), under his standing rules: commits authored as `mattlee821 <30179533+mattlee821@users.noreply.github.com>` with no `Co-Authored-By` line; before any push to a public repository, scan every tracked file and every commit for real EPIC IDs and for tokens (counts only), show Matt the file list and wait for his go-ahead.

## 1. What the project is for and what changed

The goal is to fine-map more than 10,000 GWAS (proteins from deCODE, Olink and SomaScan studies, plus others), per locus, with SuSiE-RSS and related methods, against the public Pan-UKB LD panels. The runner is `scripts/run_finemapping.py`, driven by `finemapping.nf`. The last session made the server side ready: the reference data are installed once, and the runner reads them.

## 2. Reference-data installer (done, installed on the server)

Why: the runner needs 4 reference inputs (GRCh38 to GRCh37 chain, GRCh37 FASTA, gene regions in GRCh38, and the Pan-UKB variant index of the LD panel population). They were planned but not installed. Matt chose the installer only, with its own downloads, in `~/GWASfinemapping` on the server.

What was built (all in this project, untracked):

- `scripts/reference_data.py` with `scripts/install_reference_data.sh` and `scripts/install_reference_data.sbatch`: install, check and paths for the reference root (`GWASFM_REFERENCE_ROOT`, default `~/.cache/gwasfinemapping/reference`). It installs the UCSC `hg38ToHg19.over.chain.gz`, the Ensembl release 113 GRCh37 FASTA with its own index, the gene annotation (11,105 regions from the usefuldata dataset `mapping_GRCh38_p14_uniprot`, pinned to commit `30f2e833b91c01cb0c2e22149afa89c6ba0a4ad1`; accessions that map to more than 1 gene go to `gene_regions_GRCh38.ambiguous.tsv`), and 6 whole-population Pan-UKB variant indexes (AFR, AMR, CSA, EAS, EUR, MID; GRCh37; `panukb/UKBB.<POP>.variant_index.parquet` plus a `.json` record with `region: "genome"`). `reference_manifest.json` records source, size, SHA-256 and date of every file.
- `docs/reference-data.md` documents it. `tests/python/test_reference_data.py` has 13 tests (one needs a real Rscript).
- `scripts/source.sh`: `PYTHONNOUSERSITE=1` before both pip calls and a pip upgrade, because pip 24.2 on the server failed with `ResolutionImpossible` and `~/.local` numpy and scipy hid hail's dependencies. `requirements-tools.txt` pins fsspec, s3fs, gcsfs, aiobotocore, boto3 and botocore. `Dockerfile` upgrades pip first.
- `scripts/build_panukb_variant_index.py`: Spark `driver.host` and `bindAddress` set to 127.0.0.1 (override `GWASFM_SPARK_DRIVER_HOST`), because the node hostname `cn17_mgmt` contains an underscore and Spark rejected it.

Outcome on the server (16 Nov 2026): `bash scripts/install_reference_data.sh check` ends with "all reference data are in place". Index sizes: AFR 22,082,831 variants, AMR 11,713,690, CSA 16,161,130, EAS 10,225,280, EUR 23,960,350, MID 14,127,512. Each index took about 26 min and 2.2 GB of memory (6 parallel jobs, one per population: `sbatch --job-name=gwasfm_$p scripts/install_reference_data.sbatch --only panukb --populations $p`). The reference root is 3.6 GB. Hail 0.2.135 runs under Java 21 with a warning.

## 3. Runner change (done 17 Nov 2026, tested locally)

Why: the runner built a small regional index through Hail and S3 for every locus, which is slow and needs the network, and is wasteful at 10,000 GWAS.

What changed in `scripts/run_finemapping.py`:

- New option `--panukb-dir` (default `panukb/` under `GWASFM_REFERENCE_ROOT`, else `~/.cache/gwasfinemapping/reference/panukb`).
- `ensure_regional_index` now cuts the locus rows out of the installed whole-population index with pyarrow (`cut_regional_index`, a row-group filter on contig and position) into the same regional parquet and `.json` record in `--cache-dir` that `ldcov` already read, so everything after the cut is unchanged. If the whole index is not installed it falls back to the old Hail build. An index record that is not the whole-genome GRCh37 index of that population stops the run (`installed_population_index`). A region with no variants raises an error and leaves no file.
- 4 new tests in `tests/python/test_finemapping.py` (`InstalledPopulationIndexTests`): the cut keeps only the region and the 5 columns, the cut is reused and its record is right, an empty region errors, and a missing or other index is not used. They pass. Of the 14 tests in that file, 13 pass; `test_region_method_runner_consumes_aligned_ld_matrix` fails on the Mac only because the R package jsonlite is not installed there.
- Not yet run on the server: a real locus through `ldcov` against the cut index. `ldcov` is not installed on the Mac, so the cut was tested with pyarrow only. First server check: run 1 GWAS with `--max-loci 2` and confirm the regional index is cut (a `.parquet` and `.json` in `--cache-dir/panukb/`) and that `match_variants` gives the same matches as the old regional build for the same locus.

### Server test of the runner change (to do first)

Copy the project with `bash /Users/matthewlee/research-local/copy_finemapping_to_server.sh` (Mac; it rsyncs to `~/GWASfinemapping`, excluding .git, archive, temp, .env, caches, docs/papers, parquet, data and references). Then on the server, in `~/GWASfinemapping` with the tools environment active (`bash scripts/source.sh --tools-only` has been run, and `bash scripts/install_reference_data.sh check` ends with "all reference data are in place"), run 1 GWAS-VCF (for example 1 of the 18 Sun files in `drive1/research/GWASprocessing/data/GWAS/sun_2023_37794186/standardise/`, copied to the server) with `--max-loci 2` and the paths from `bash scripts/install_reference_data.sh paths --nextflow`. Confirm 3 things: (1) `--cache-dir/panukb/UKBB.<POP>.chr<chrom>_<start>_<end>.variant_index.parquet` and its `.json` appear for each locus, with no Hail or Spark output in the log; (2) the number of variants matched to the panel (`match_variants`) equals the number from a run with `--panukb-dir` pointing at an empty folder, which forces the old regional Hail build; (3) the run time of the index step is seconds. A difference in (2) is a bug in `cut_regional_index`. Peak memory of `ldcov` reading a cut index should be checked on the largest locus (MHC).

## 4. Pieces of the pipeline outside this project that this one uses

- GWAS inputs are GWAS-VCFs from `~/GWASprocessing` (standardise step). For the first fine-mapping runs, the 18 Sun UKB-PPP Olink GWAS (12 proteins) are on drive1 under `drive1/research/GWASprocessing/data/GWAS/sun_2023_37794186/standardise/`.
- A cis-variant clump for those 12 proteins was run from `temp/` (untracked, own `.gitignore`): `cis_clump.R`, `cis_clump.sbatch`, `build_ld_reference.sh`. Gene region plus 1 Mb either side, p < 5e-8, plink2 `--clump` at r2 0.1, 0.01 and 0.001 against the 1000 Genomes GRCh38 30x phase 3 unrelated samples (EUR, 503; ALL, 2,504). Lessons: the 1000 Genomes rsIDs repeat within a chromosome, so IDs are set to `chr:pos:ref:alt` (long alleles truncated, exact repeats dropped); 63% of the EUR reference variants were monomorphic and counted as independent, so the reference keeps only MAF at least 0.01 (a GWAS variant below that has NA in the flag columns); r2 0.001 on 503 samples is indicative only. IL6 (4 GWAS) and IFNG have no cis variant at p < 5e-8. Tables and a README are on drive1 at `.../sun_2023_37794186/cis_clump/`.

## 5. Step to redo (the agent was mid-way through it)

The previous agent working on this project had just said: "I’ll review the new paper and current output/locus design, then recommend a storage-light output contract and defensible region/MHC defaults before changing anything." That was in answer to the following message from Matt. Redo this step from the start, and make no code change until Matt has agreed the recommendation.

Matt's message (verbatim):

"i dont think we need to keep all of these outputs; remember that we will have >10,000 finemapped GWAS so we need to be concious of space. would it be possible to only keep the .rda results files and the .json files - this would give us all the information we required for downstream analysis wouldnt it? this might require us to create a report.html which summarises all of the results files including flags and displays figures for each output. what do you think?
we need to make a decision on default values for clusters; i cant see any information in the literature review for clustering or the window size; i have added a new paper which defined regions using a 3Mb window; we also need to make a decision about the MHC region, i think this should be included as a single region by default and flags should be added (e.g., so downstream exclude MHC can be applied) but we should also be able to exclude the MHC when running the pipeline. can we think about how we can use LD defined boundaries for the finemapping window? or perhaps it would be best if we copied the GeneHackman approach if there is no clear answer from the literature review"

What the step has to produce, as a recommendation for Matt to approve:

1. Storage-light output contract. Decide whether the `.rda` results files and the `.json` files hold everything downstream analysis needs (credible sets, PIPs, per-variant statistics, flags, LD and method settings, run provenance), name what they lack, and list every other output the runner writes now (per locus: lifted VCF, normalised variants, mapping tables, LD matrices and npz, method folders, png figures, tsv) with what can be deleted after a run. Estimate bytes per GWAS and per 10,000 GWAS before and after. Design a `report.html` per GWAS (or per batch) built from the `.rda` and `.json` files: summary of loci, flags and 1 figure set per locus, with no need to keep the figure files.
2. Locus clustering and window defaults. The current defaults in `scripts/run_finemapping.py` are `--significance-lp` at p 5e-8, `--merge-distance-bp` 250,000, `--region-flank-bp` 100,000 and `--gene-window-bp` 1,000,000 (gene-region analysis). `docs/literature_review.md` (item 1 of its open decisions, and the GeneHackman reference at `R/finemap.R`) names these as needing evidence-backed defaults but gives none. Matt added a paper that defined regions with a 3 Mb window: it is `docs/papers/2021.09.03.21262975v1.full.pdf` (confirmed by Matt on 17 Nov 2026). Read its full text, report exactly how it defined and merged regions, and say which parts transfer to protein GWAS (cis and trans) and to colon and rectal cancer GWAS. Do not cite from an abstract or search summary, and say which papers could not be obtained.
3. MHC region (chr6, GRCh37 about 25–34 Mb, GRCh38 about 28.5–33.4 Mb; confirm the coordinates from a source). Matt's position: include the MHC as a single region by default, add a flag on every MHC-overlapping locus and result so downstream analyses can exclude it, and add a pipeline option that excludes the MHC from the run. Check this against the literature and the LD panel (Pan-UKB LD in the MHC) and say what the single-region window does to runtime and memory.
4. LD-defined boundaries for the window. Assess whether LD-block boundaries (for example published LD blocks, or boundaries from the Pan-UKB LD panel) are a defensible window for the fine-mapping region, compare with a fixed window (1 Mb or 3 Mb) and with the GeneHackman approach, and recommend one. If the literature gives no clear answer, copy GeneHackman as Matt suggested and say so.
5. Where the decisions are written: reusable methods decisions go in the methods skill (`~/.claude/skills/methods/references/`, one reference file per topic); choices specific to this project go in this project's plan or `docs/`.

## 6. Open items for this project

- Test the runner change on the server (section 3) before any large run.
- The 5 Pan-UKB populations other than EUR are installed but the 5 other LD panels have not been tested with a real locus.
- The nextflow wrapper `finemapping.nf` takes the reference paths from `bash scripts/install_reference_data.sh paths --nextflow`; check that it passes `--panukb-dir`, or relies on the default under `GWASFM_REFERENCE_ROOT`.
- Committing and pushing the restructured project is the job of the agent on this project (see the top of this note for the rules).
