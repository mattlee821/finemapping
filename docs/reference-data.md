# Reference data

A fine-mapping run needs 4 reference inputs besides the GWAS: the GRCh38 to GRCh37 chain and the GRCh37 FASTA (to move
the GWAS variants onto the Pan-UKB LD panel), the gene annotation (gene regions in GRCh38), and the Pan-UKB variant index
of the LD panel population. `scripts/install_reference_data.sh` installs them once, under 1 reference root.

## What is installed

The reference root is `GWASFM_REFERENCE_ROOT` in `.env` (default `~/.cache/gwasfinemapping/reference`, outside
cloud-synced folders).

| Path under the root | Content | Source |
|---|---|---|
| `genome/hg38ToHg19.over.chain.gz` | GRCh38 to GRCh37 chain | UCSC (`hg38ToHg19.over.chain.gz`) |
| `genome/GRCh37.fa`, `GRCh37.fa.fai` | GRCh37 primary assembly (about 3 GB), contigs `1`-`22`, `X`, `Y` | Ensembl release 113, indexed by the installer |
| `annotation/gene_regions_GRCh38.tsv` | 1 gene region per UniProt accession (the layout of `config/gene_regions_test.tsv`, plus `source_rows`) | the dataset `mapping_GRCh38_p14_uniprot` of the public R package [usefuldata](https://github.com/mattlee821/usefuldata), pinned to a commit |
| `annotation/gene_regions_GRCh38.ambiguous.tsv` | accessions whose rows name more than 1 gene, left out of the annotation (the runner needs exactly 1 gene per accession) | same |
| `polyfun/by_chrom/chr1..22.parquet` and `polyfun_meta.json` | 19,476,620 autosomal GRCh37 per-SNP functional heritability priors (`baselineLF2.2.UKB` meta-analysis across 15 UK Biobank traits) with genome-wide $q = 100$ max/min truncation (`4.0733e-09` to `4.0733e-07`) | [omerwe/polyfun](https://github.com/omerwe/polyfun) (Weissbrod et al. 2020), pinned to commit `8b17b4ec` |
| `panukb/UKBB.<POP>.variant_index.parquet` and `.json` | the whole-population variant index of the public Pan-UKB LD panel (GRCh37), for AFR, AMR, CSA, EAS, EUR and MID; `.json` records the row count, source table and software versions | `s3://pan-ukb-us-east-1/ld_release/UKBB.<POP>.ldadj.variant.ht`, read with Hail |
| `reference_manifest.json` | source, size, SHA-256 and date of every installed file | written by the installer |

Rows of the same gene in the usefuldata mapping are merged into 1 span (minimum start, maximum end), so a gene with
several annotated copies gets 1 wide region. Alternative contigs and rows without a GRCh38 position are left out.

## Install on a server

```bash
# on the server, in the project folder (for example ~/GWASfinemapping), after the project files are there
cp -n .env.example .env                    # then set GWASFM_REFERENCE_ROOT if the default is not where you want it
bash scripts/source.sh --tools-only        # micromamba environment: Python, Java, Hail, ldcov, R (needs network)
bash scripts/install_reference_data.sh --only genome,annotation   # small and quick: run it on the login node
sbatch scripts/install_reference_data.sbatch                      # everything, including the 6 Pan-UKB indexes
```

The job installs everything and then runs `check`. Each step skips what is already installed and valid, so a failed job
can be submitted again. Options (also for `bash scripts/install_reference_data.sh install ...`):

- `--only genome,annotation,panukb` limits the steps; `--populations EUR,AFR` limits the Pan-UKB indexes.
- `--force` reinstalls what is already there; `--dry-run` shows the plan.
- `GWASFM_CHAIN_URL`, `GWASFM_FASTA_GRCH37_URL` and `GWASFM_USEFULDATA_COMMIT` change the sources.
- `GWASFM_HAIL_DRIVER_MEMORY` (default `8g` here) and `GWASFM_HAIL_MASTER` (default `local[4]`) size the Hail session that
  reads the Pan-UKB table.

Each Pan-UKB index takes about 26 min and 2.2 GB of memory, so the 6 can run as 6 parallel jobs (`sbatch --job-name=gwasfm_$p scripts/install_reference_data.sbatch --only panukb --populations $p`). Installed on the server on 16 Nov 2026: 3.6 GB in all, indexes of 10.2 million (EAS) to 24.0 million (EUR) variants.

The compute node needs network access to UCSC, Ensembl, GitHub and the public S3 bucket.

## Check and use

```bash
bash scripts/install_reference_data.sh check                 # what is missing or changed; exit status 1 if anything is
bash scripts/install_reference_data.sh paths --nextflow      # --chain ... --grch37_fasta ... --grch37_fai ... --gene_annotation ...
nextflow run finemapping.nf $(bash scripts/install_reference_data.sh paths --nextflow) --vcf ... --population EUR ...
```

`check` re-reads each file: the chain must run from GRCh38 chromosome 1 (248,956,422 bases) to GRCh37 (249,250,621), the
FASTA must hold GRCh37 contigs and be complete, the annotation must have 1 valid GRCh38 region per accession, and each
Pan-UKB index must have the row count its record states.

## How the runner uses it

`scripts/run_finemapping.py` finds the whole-population index of its `--population` in `--panukb-dir` (default: `panukb/` under
`GWASFM_REFERENCE_ROOT`, else `~/.cache/gwasfinemapping/reference/panukb`) and cuts the rows of each locus out of it with
pyarrow into a small regional index in `--cache-dir`, which `ldcov` then reads as before. No Hail session and no network are
needed, and the cut takes seconds. If the index is not installed, the runner builds the regional index through Hail and S3
as it did before. An index whose record is not the whole-genome GRCh37 index of the population stops the run.
