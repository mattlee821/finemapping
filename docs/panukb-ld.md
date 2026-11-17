# Pan-UKB LD from public AWS S3

The Pan-UKB LD BlockMatrices and companion variant-index Hail Tables are publicly
readable from AWS S3. The Google Cloud Storage convenience bundle used by `ldcov`
is requester-pays, but it is not required for this route.

The released LD matrices are ancestry-specific: `AFR`, `AMR`, `CSA`, `EAS`, `EUR`,
and `MID`. There is no `ALL` BlockMatrix. A pooled or multi-ancestry GWAS must not
be treated as ancestry-matched to one of these panels without a justified analysis
method.

## Install the tools

The project setup follows GWASprocessing's micromamba pattern. It installs a
private environment outside cloud-synced folders by default, with Python 3.12,
Java 21, Hail, `ldcov[s3]`, CMake, and the R package dependencies from
`renv.lock`.
On Intel Macs, setup builds `lazybgen` with its non-assembly zstd path because
the upstream package does not publish an Intel macOS wheel and its default source
build fails to link there.

```bash
bash scripts/source.sh
source scripts/source.sh
```

Use `bash scripts/source.sh --tools-only` to omit the R package restore. The setup
pins Hail and `ldcov` in [requirements-tools.txt](../requirements-tools.txt). Its
optional container mode pulls the published GHCR image with `--container docker`
or `--container apptainer`.

## Build a local Parquet variant index

The project builder configures Hail/Spark's Hadoop S3A connector for anonymous
access to the public bucket. It downloads the Hadoop AWS connector and dependencies
from Maven Central on first use.

```bash
source scripts/source.sh
bash scripts/build_panukb_variant_index.sh \
  EUR \
  /path/to/UKBB.EUR.chr22_25000000_25100000.variant_index.parquet \
  22:25000000-25100000
```

The script selects the public S3 table
`s3a://pan-ukb-us-east-1/ld_release/UKBB.EUR.ldadj.variant.ht` and exports the
columns required by `ldcov`. Choose one of the six released populations explicitly.
Set `GWASFM_HAIL_MASTER` and `GWASFM_HAIL_DRIVER_MEMORY` to adjust the local Spark
master and driver memory. The optional region argument limits the index to that
GRCh37 region; the resulting Parquet index should only be used for extracting LD
for that same region. Omit it to create a whole-population index.

The builder writes the selected Hail Table rows to a temporary compressed TSV, then
streams them into Parquet in 250,000-row chunks. It avoids collecting all variants
in Pandas at once. A whole-population index still scans and temporarily exports the
complete variant table, but does not download the LD BlockMatrix.

## Extract a regional LD matrix

Request a small region from the public BlockMatrix. The input index and region
coordinates must both be GRCh37:

```bash
bash scripts/extract_panukb_ld_region.sh \
  EUR \
  /path/to/UKBB.EUR.chr22_25000000_25100000.variant_index.parquet \
  22:25000000-25100000 \
  /path/to/results/eur_chr22_test
```

This writes an `.npz` matrix and a `.variants.tsv` file; it does not download the
full BlockMatrix. Public access avoids the GCS requester-pays bundle, but local
compute, storage, and network-transfer constraints still apply.

The project GWAS-VCFs are GRCh38. These scripts intentionally do not perform
GRCh38-to-GRCh37 liftover or allele matching. Do not pass GRCh38 coordinates
directly: an allele-aware build mapping and variant-order validation are still
required before using extracted LD for fine-mapping.

## Locke BMI ancestry

The Locke input metadata labels the GWAS ancestry `ALL` and identifies the source
as an all-ancestry summary-statistics file. Its sample size also varies by variant.
Using UK Biobank's ancestry-specific LD panel for those pooled statistics would not
create a matching combined-ancestry LD reference. Keep Locke to P-value
prioritization until ancestry-stratified summary statistics and matching LD are
available, or until a suitable multi-ancestry fine-mapping method and all required
inputs are selected. Do not infer European ancestry from the 1000 Genomes EUR
frequency annotation used to fill missing allele frequencies.

## Bounded Nextflow directory preflight pilot

The pilot discovers standardized GWAS-VCF directories, writes a manifest, and
runs one lightweight preflight task per study. It checks metadata/QC build
agreement, the VCF header and sample column, and the tabix index record count
using `bcftools index --nrecords`. It does not scan all VCF records, copy the
large indexed inputs into task work directories, perform liftover, retrieve LD,
or run fine-mapping. `bcftools` is installed in the project tools environment.

Run the pilot with Nextflow 23.10 or later:

```bash
source scripts/source.sh
nextflow run main.nf \
  --input_dir /path/to/ferkingstad_2021_34857953/standardise \
  --outdir "$HOME/.cache/gwasfinemapping/ferkingstad-preflight" \
  --max_inputs 10 \
  --max_forks 2
```

`--max_inputs` defaults to 10; set it to `0` only to select every discovered
directory. `--max_forks` bounds concurrent study tasks. The default local
executor can be changed with `-profile slurm`; the input and work paths must then
be visible on the compute nodes. The work directory defaults to
`~/.cache/gwasfinemapping/nextflow-work`, outside cloud-synced folders.

The output contains the selected `manifest.tsv`, `discovery_summary.json`, one
JSON report per study, and `pilot_summary.tsv`. A new or empty output directory
is required for a fresh run. To resume, reuse that output directory with
Nextflow's `-resume`; the workflow rechecks that the selected manifest is
unchanged and refuses a mismatched run. The current local pilot completed on the
two available Ferkingstad directories: their tabix indexes report 25,359,358
and 25,359,359 records, respectively. Their compressed VCFs are about 849 MB
each. This
validates manifest handling, symlink staging, and bounded task execution; it is
not a 10–100-study throughput benchmark. Only two directories were present at
the specified location when tested.

## Bounded end-to-end fine-mapping test

The one-GWAS Nextflow workflow completed on the Ferkingstad P43320/CRYBB2
GWAS-VCF using the EUR public-S3 LD panel. A one-locus run and a bounded
three-locus run both completed. The three-locus run matched 893, 897, and 82
variants to LD; each resulting matrix was complete and symmetric, with raw
diagonal values close to one. P-value, p-LD, SuSiE-RSS, and SuSiE-rRSS outputs
were produced for each locus. The CRYBB2 gene-body region extended by 1 Mb
contained 20,280 variants, all labelled cis.

## Result layout

- `discovery/` contains the genome-wide significant-variant table and the
  clustered-locus table, including lead-variant liftover status.
- `gene_region/` contains the gene-body ± cis-window variants and their
  prioritization sidecars.
- `loci/<region-id>/` contains per-locus inputs, variant-mapping audit,
  validated LD, and method-specific TSV/JSON sidecars.
- `methods/` contains one aggregate R workspace per method: `p-value.rda`,
  `p-ld.rda`, `susie-rss.rda`, and `susie-rrss.rda`. Each file contains one
  object named for its method ID (for example, `susie_rss`), with a `regions`
  list. Each region entry has `region_id`, `region_type`, and `result`; the
  result retains method outputs, diagnostics, provenance, and (for SuSiE) the
  fitted model. The gene-region and genome-wide-locus results are distinguished
  by `region_type`.

For example, load the aggregate SuSiE-RSS output with
`load("methods/susie-rss.rda")`, list the analyzed regions with
`names(susie_rss$regions)`, then inspect a selected
`susie_rss$regions[[region_id]]$result`. TSV and JSON
sidecars remain available for quick inspection and audit without R.

## Locus selection and fine-mapping intervals

The runner selects genome-wide significant variants at `LP >= 7.30103`
(`P <= 5e-8` by default), sorts them by chromosome and position, and clusters
adjacent significant variants when their positional gap is at most 250 kb.
This is adjacency chaining: a cluster can span more than 250 kb when successive
significant variants bridge the gaps. The cluster lead is its most significant
variant. The lead must map to GRCh37 for the cluster to be eligible for
fine-mapping.

The runner now analyzes every mapped cluster by default. Set `--max-loci N` to
cap a smoke test; `--max-loci 0` means no cap. For each eligible cluster, the
GWAS query and LD matching interval spans the full cluster start/end plus the
configured flank on both sides (100 kb by default). This prevents the analysis
window from excluding significant cluster members far from the lead. The
250-kb clustering distance and 100-kb flank remain provisional settings, not
validated biological or LD boundaries.

These rules differ from GeneHackman's workflow. The documented example uses
PLINK clumping at `P < 5e-8`, `r² = 0.01`, and a 1-Mb clumping distance, then
fine-maps lead-centered windows of about ±1 Mb and merges overlapping windows
(up to a 3-Mb region). This project groups by positional distance without LD
clumping, and then analyzes the cluster span plus its flank. Therefore its locus
list and intervals are not directly comparable to GeneHackman's. The prior
controlled comparison held variants and LD fixed to compare the SuSiE fitting
step only; it did not compare the complete workflows.

This verifies the end-to-end data path, not the scientific reliability of every
fine-mapping result. At the strongest mapped locus, SuSiE-RSS did not converge
with either 100 or 200 iterations. SuSiE-rRSS did not converge at 100, but
converged at iteration 124 with a 200-iteration cap; its LD/summary-statistics
reliability warning remained. `SER_fallback` is set equal to the reliability
flag; it does not replace the multi-effect fit. A one-effect diagnostic is
available within `R_finite_diagnostics` for inspection. These diagnostics are
retained with the outputs and must be investigated before interpreting that
locus. That historical test explicitly capped analysis at three loci and used
lead-centered windows; the current runner uses cluster-span windows and defaults
to all mapped clusters. The significance threshold, positional clustering
distance, interval flank, and cis window remain provisional settings, not
validated production choices. The workflow accepts one GWAS per run;
directory-scale scheduling and throughput have not been tested.

## Validation status

Public S3 access was verified with Hail and `ldcov`. Building the EUR index for
GRCh37 `22:25000000-25100000` returned 700 variants; extraction returned a
700-by-700 symmetric LD matrix whose diagonal ranged from 0.999919 to 0.999989.
This validates the public-S3 regional extraction path. The bounded end-to-end
test above additionally validates allele-aware GRCh38-to-GRCh37 matching and LD
consumption for one GWAS, subject to the reported method diagnostics. The
full-population index build and directory-scale fine-mapping remain untested.
