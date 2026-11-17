#!/usr/bin/env python3
"""Install and check the reference data that GWASfinemapping needs at run time.

Usage:
    python3 scripts/reference_data.py install [--root DIR] [--only genome,panukb,annotation]
                                              [--populations EUR,AFR,AMR,CSA,EAS,MID] [--force] [--dry-run]
    python3 scripts/reference_data.py check   [--root DIR] [--only ...] [--populations ...]
    python3 scripts/reference_data.py paths   [--root DIR] [--nextflow]

The reference root (--root, else GWASFM_REFERENCE_ROOT, else ~/.cache/gwasfinemapping/reference, which is outside
cloud-synced folders) holds:

    genome/hg38ToHg19.over.chain.gz      UCSC GRCh38 to GRCh37 chain, for the liftover of the GWAS variants
    genome/GRCh37.fa, GRCh37.fa.fai      Ensembl GRCh37 primary assembly, to normalise the lifted alleles
    panukb/UKBB.<POP>.variant_index.parquet  (+ .json)
                                         the whole-population variant index of each public Pan-UKB LD panel (GRCh37),
                                         read from the public S3 bucket with Hail; needs the tools environment
    annotation/gene_regions_GRCh38.tsv   1 gene region per UniProt accession from usefuldata (GRCh38), in the layout of
                                         config/gene_regions_test.tsv; accessions that map to more than 1 gene are written
                                         to gene_regions_GRCh38.ambiguous.tsv and left out
    reference_manifest.json              where each file came from, its size and SHA-256, and when

Every step skips what is already installed and valid; --force reinstalls. `check` re-reads the files and reports what is
missing or does not match its record. Override the sources with GWASFM_CHAIN_URL, GWASFM_FASTA_GRCH37_URL and
GWASFM_USEFULDATA_COMMIT.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path
from typing import Any

POPULATIONS = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")
STEPS = ("genome", "panukb", "annotation")
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
FASTA_GRCH37_URL = ("https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/"
                    "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz")
USEFULDATA_COMMIT = "30f2e833b91c01cb0c2e22149afa89c6ba0a4ad1"  # mattlee821/usefuldata, includes the UniProt mapping
USEFULDATA_URL = "https://raw.githubusercontent.com/mattlee821/usefuldata/{commit}/data/mapping_GRCh38_p14_uniprot.rda"
USEFULDATA_DATASET = "mapping_GRCh38_p14_uniprot"
GRCH37_CHR1_LENGTH = 249_250_621
GRCH38_CHR1_LENGTH = 248_956_422
STANDARD_CONTIGS = {str(c) for c in range(1, 23)} | {"X", "Y"}
ANNOTATION_COLUMNS = ["uniprot_id", "gene_symbol", "ensembl_id", "chrom", "start", "end", "genome_build", "source", "version",
                      "source_rows"]
SCRIPTS_DIR = Path(__file__).resolve().parent


class ReferenceError(RuntimeError):
    """An installed or downloaded file is missing or not what it should be."""


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def reference_root(value: str | None = None) -> Path:
    root = value or os.environ.get("GWASFM_REFERENCE_ROOT") or str(Path.home() / ".cache" / "gwasfinemapping" / "reference")
    return Path(root).expanduser().resolve()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# manifest ====
def manifest_path(root: Path) -> Path:
    return root / "reference_manifest.json"


def read_manifest(root: Path) -> dict[str, Any]:
    path = manifest_path(root)
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    return {"schema_version": 1, "items": {}}


def record(root: Path, relative: str, **fields: Any) -> None:
    manifest = read_manifest(root)
    path = root / relative
    manifest["items"][relative] = {"bytes": path.stat().st_size, "sha256": sha256_file(path), "installed": now(), **fields}
    manifest["updated"] = now()
    root.mkdir(parents=True, exist_ok=True)
    temporary = manifest_path(root).with_suffix(".json.part")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, manifest_path(root))


# download ====
def fetch(url: str, destination: Path, retries: int = 3, timeout: int = 120) -> Path:
    """Download to <destination>.part and rename when complete; retried with a growing pause."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_name(destination.name + ".part")
    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response, part.open("wb") as out:
                shutil.copyfileobj(response, out, 8 * 1024 * 1024)
            if part.stat().st_size == 0:
                raise ReferenceError(f"{url} returned an empty file")
            os.replace(part, destination)
            return destination
        except Exception as error:  # noqa: BLE001 - retried, then raised with the cause
            last = error
            print(f"  download attempt {attempt}/{retries} failed: {error}", file=sys.stderr)
            time.sleep(min(30, 5 * attempt))
    part.unlink(missing_ok=True)
    raise ReferenceError(f"could not download {url}: {last}")


# genome ====
def build_fai(fasta: Path) -> int:
    """Write <fasta>.fai (name, length, offset, bases per line, bytes per line) like `samtools faidx`; the lines of a
    record must all have the same width except the last. Returns the number of records."""
    entries: list[tuple[str, int, int, int, int]] = []
    name = None
    length = offset = linebases = linewidth = 0
    short_seen = False
    position = 0
    with fasta.open("rb") as stream:
        for raw in stream:
            width = len(raw)
            text = raw.rstrip(b"\r\n")
            if raw.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, linebases, linewidth))
                name = text[1:].split()[0].decode()
                length, offset, linebases, linewidth, short_seen = 0, position + width, 0, 0, False
            elif name is None:
                raise ReferenceError(f"{fasta}: sequence before the first header")
            elif text:
                if linebases == 0:
                    linebases, linewidth = len(text), width
                elif short_seen or len(text) > linebases:  # only the last line of a record may be shorter
                    raise ReferenceError(f"{fasta}: uneven line widths in {name}")
                if len(text) < linebases:
                    short_seen = True
                length += len(text)
            position += width
    if name is not None:
        entries.append((name, length, offset, linebases, linewidth))
    if not entries:
        raise ReferenceError(f"{fasta}: no records")
    with Path(f"{fasta}.fai").open("w", encoding="utf-8") as out:
        for entry in entries:
            out.write("\t".join(str(v) for v in entry) + "\n")
    return len(entries)


def check_fasta(fasta: Path) -> None:
    fai = Path(f"{fasta}.fai")
    if not fasta.is_file() or not fai.is_file():
        raise ReferenceError(f"{fasta} or its .fai is missing")
    lengths = {row[0]: int(row[1]) for row in (line.split("\t") for line in fai.read_text().splitlines())}
    missing = sorted(STANDARD_CONTIGS - set(lengths))
    if missing:
        raise ReferenceError(f"{fasta} lacks contigs {', '.join(missing)} (names must be 1-22, X, Y without 'chr')")
    if lengths["1"] != GRCH37_CHR1_LENGTH:
        raise ReferenceError(f"{fasta}: chromosome 1 has {lengths['1']} bases, not {GRCH37_CHR1_LENGTH} (not GRCh37)")
    if sum(lengths.values()) < 2_800_000_000:
        raise ReferenceError(f"{fasta}: {sum(lengths.values())} bases in all; the file looks truncated")


def check_chain(chain: Path) -> None:
    """The UCSC hg38ToHg19 chain: a gzip file of chain records from GRCh38 (target) to GRCh37 (query)."""
    if not chain.is_file():
        raise ReferenceError(f"{chain} is missing")
    found = False
    try:
        with gzip.open(chain, "rt", encoding="utf-8") as stream:
            for line in stream:
                if line.startswith("chain "):
                    fields = line.split()
                    if fields[2] == "chr1" and int(fields[3]) == GRCH38_CHR1_LENGTH and fields[7] == "chr1" and int(fields[8]) == GRCH37_CHR1_LENGTH:
                        found = True
                        break
    except (OSError, EOFError, ValueError, IndexError) as error:
        raise ReferenceError(f"{chain} cannot be read as a chain file: {error}") from error
    if not found:
        raise ReferenceError(f"{chain} has no chr1 chain from GRCh38 ({GRCH38_CHR1_LENGTH}) to GRCh37 ({GRCH37_CHR1_LENGTH})")


def install_genome(root: Path, force: bool) -> None:
    out = root / "genome"
    chain = out / "hg38ToHg19.over.chain.gz"
    url = os.environ.get("GWASFM_CHAIN_URL", CHAIN_URL)
    if force or not chain.is_file():
        print(f"genome: downloading the chain from {url}")
        fetch(url, chain)
    check_chain(chain)
    record(root, "genome/hg38ToHg19.over.chain.gz", source=url)

    fasta = out / "GRCh37.fa"
    fasta_url = os.environ.get("GWASFM_FASTA_GRCH37_URL", FASTA_GRCH37_URL)
    if force or not fasta.is_file():
        print(f"genome: downloading the GRCh37 FASTA from {fasta_url} (about 0.9 GB, 3 GB unpacked)")
        compressed = out / "GRCh37.fa.gz"
        fetch(fasta_url, compressed)
        part = out / "GRCh37.fa.part"
        with gzip.open(compressed, "rb") as source, part.open("wb") as target:
            shutil.copyfileobj(source, target, 8 * 1024 * 1024)
        os.replace(part, fasta)
        compressed.unlink()
        Path(f"{fasta}.fai").unlink(missing_ok=True)
    if not Path(f"{fasta}.fai").is_file():
        print(f"genome: indexed {build_fai(fasta)} records")
    check_fasta(fasta)
    record(root, "genome/GRCh37.fa", source=fasta_url)
    record(root, "genome/GRCh37.fa.fai", source="built from GRCh37.fa")


# Pan-UKB variant indexes ====
def panukb_paths(root: Path, population: str) -> tuple[Path, Path]:
    index = root / "panukb" / f"UKBB.{population}.variant_index.parquet"
    return index, Path(f"{index}.json")


def check_panukb(root: Path, population: str) -> int:
    index, sidecar = panukb_paths(root, population)
    if not index.is_file() or not sidecar.is_file():
        raise ReferenceError(f"{index} or its record is missing")
    expected = json.loads(sidecar.read_text(encoding="utf-8"))
    if (expected.get("population"), expected.get("region"), expected.get("build")) != (population, "genome", "GRCh37"):
        raise ReferenceError(f"{sidecar} describes another index")
    import pyarrow.parquet as pq

    parquet = pq.ParquetFile(index)
    rows = parquet.metadata.num_rows
    missing = {"contig", "position", "ref", "alt", "idx"} - set(parquet.schema_arrow.names)
    if rows == 0 or missing or rows != expected.get("rows"):
        raise ReferenceError(f"{index}: {rows} rows (recorded {expected.get('rows')}), missing columns {sorted(missing)}")
    return rows


def install_panukb(root: Path, populations: tuple[str, ...], force: bool) -> None:
    for population in populations:
        index, sidecar = panukb_paths(root, population)
        if not force:
            try:
                rows = check_panukb(root, population)
                print(f"panukb {population}: installed ({rows:,} variants)")
                continue
            except ReferenceError:
                pass
        for stale in (index, sidecar, Path(f"{index}.tmp"), Path(f"{index}.export.tsv.bgz")):
            if stale.exists():  # the builder refuses to overwrite; an earlier interrupted build leaves these
                print(f"panukb {population}: removing {stale.name}")
                stale.unlink()
        env = {**os.environ, "GWASFM_HAIL_DRIVER_MEMORY": os.environ.get("GWASFM_HAIL_DRIVER_MEMORY", "8g"),
               "GWASFM_HAIL_MASTER": os.environ.get("GWASFM_HAIL_MASTER", "local[4]")}
        print(f"panukb {population}: building the whole-population index from the public S3 table (this takes a while)")
        started = time.time()
        subprocess.run(["bash", str(SCRIPTS_DIR / "build_panukb_variant_index.sh"), population, str(index)], check=True, env=env)
        rows = None
        try:
            import pyarrow.parquet as pq

            rows = pq.ParquetFile(index).metadata.num_rows
        except Exception as error:  # noqa: BLE001
            raise ReferenceError(f"{index} was built but cannot be read: {error}") from error
        versions = {}
        try:
            from importlib.metadata import version

            versions = {name: version(name) for name in ("hail", "ldcov", "lazybgen")}
        except Exception:  # noqa: BLE001 - versions are provenance only
            pass
        sidecar.write_text(json.dumps({
            "population": population, "region": "genome", "build": "GRCh37", "rows": rows,
            "source": f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.variant.ht", "built": now(),
            "seconds": round(time.time() - started), "software": versions}, indent=2) + "\n", encoding="utf-8")
        check_panukb(root, population)
        record(root, f"panukb/{index.name}", source=f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.variant.ht", rows=rows)
        print(f"panukb {population}: {rows:,} variants")


# gene annotation ====
def build_gene_annotation(table: list[dict[str, str]], version: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """1 region per UniProt accession from the usefuldata mapping rows (UNIPROT, uniprot_gn_symbol, ensembl_gene_id, CHR,
    START_hg38, END_hg38). Rows of the same gene (same symbol and chromosome) are merged into 1 span, with the first
    Ensembl ID there is. An accession whose rows name more than 1 gene is ambiguous and returned separately."""
    by_accession: dict[str, dict[tuple[str, str], list[dict[str, str]]]] = {}
    for row in table:
        accession, symbol = (row.get("UNIPROT") or "").strip(), (row.get("uniprot_gn_symbol") or "").strip()
        chrom, start, end = (row.get("CHR") or "").strip(), (row.get("START_hg38") or "").strip(), (row.get("END_hg38") or "").strip()
        if not accession or not symbol or chrom not in STANDARD_CONTIGS or not start.isdigit() or not end.isdigit():
            continue
        by_accession.setdefault(accession, {}).setdefault((symbol, chrom), []).append(row)
    annotation, ambiguous = [], []
    for accession in sorted(by_accession):
        genes = by_accession[accession]
        if len(genes) > 1:
            ambiguous.append({"uniprot_id": accession, "genes": ";".join(f"{s}@{c}" for s, c in sorted(genes))})
            continue
        (symbol, chrom), rows = next(iter(genes.items()))
        ensembl = next((r["ensembl_gene_id"] for r in rows if (r.get("ensembl_gene_id") or "").strip()), "")
        annotation.append({
            "uniprot_id": accession, "gene_symbol": symbol, "ensembl_id": ensembl, "chrom": chrom,
            "start": min(int(r["START_hg38"]) for r in rows), "end": max(int(r["END_hg38"]) for r in rows),
            "genome_build": "GRCh38", "source": "usefuldata", "version": version, "source_rows": len(rows)})
    return annotation, ambiguous


def export_usefuldata(rda: Path, destination: Path, rscript: str) -> None:
    code = ("args <- commandArgs(TRUE); e <- new.env(); load(args[1], e); x <- get(args[3], e); "
            "write.table(x, args[2], sep = '\\t', quote = FALSE, row.names = FALSE, na = '')")
    result = subprocess.run([rscript, "--vanilla", "-e", code, str(rda), str(destination), USEFULDATA_DATASET],  # Rscript adds --args itself
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise ReferenceError(f"Rscript could not read {rda}: {result.stderr.strip()[:300]}")


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    with temporary.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def install_annotation(root: Path, force: bool, rda: Path | None = None, rscript: str | None = None) -> None:
    out = root / "annotation" / "gene_regions_GRCh38.tsv"
    ambiguous_out = root / "annotation" / "gene_regions_GRCh38.ambiguous.tsv"
    if out.is_file() and not force:
        check_annotation(root)
        print(f"annotation: installed ({out})")
        return
    commit = os.environ.get("GWASFM_USEFULDATA_COMMIT", USEFULDATA_COMMIT)
    with tempfile.TemporaryDirectory(prefix="gwasfm_usefuldata_") as work:
        work_path = Path(work)
        if rda is None:
            url = USEFULDATA_URL.format(commit=commit)
            print(f"annotation: downloading {USEFULDATA_DATASET} from usefuldata {commit[:8]}")
            rda = fetch(url, work_path / f"{USEFULDATA_DATASET}.rda")
        exported = work_path / "mapping.tsv"
        export_usefuldata(rda, exported, rscript or shutil.which("Rscript") or "Rscript")
        with exported.open(encoding="utf-8", newline="") as stream:
            table = list(csv.DictReader(stream, delimiter="\t"))
    annotation, ambiguous = build_gene_annotation(table, f"{USEFULDATA_DATASET} {commit}")
    if len(annotation) < 5000:
        raise ReferenceError(f"only {len(annotation)} gene regions came out of {len(table)} mapping rows; the mapping looks wrong")
    write_tsv(out, annotation, ANNOTATION_COLUMNS)
    write_tsv(ambiguous_out, ambiguous, ["uniprot_id", "genes"])
    record(root, "annotation/gene_regions_GRCh38.tsv", source=f"usefuldata {commit}", regions=len(annotation), ambiguous=len(ambiguous))
    print(f"annotation: {len(annotation):,} gene regions; {len(ambiguous):,} accessions map to more than 1 gene and are left out")


def check_annotation(root: Path) -> int:
    path = root / "annotation" / "gene_regions_GRCh38.tsv"
    if not path.is_file():
        raise ReferenceError(f"{path} is missing")
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    seen = set()
    for row in rows:
        if row["genome_build"] != "GRCh38" or int(row["start"]) < 1 or int(row["end"]) < int(row["start"]) or row["uniprot_id"] in seen:
            raise ReferenceError(f"{path}: invalid or repeated region for {row['uniprot_id']}")
        seen.add(row["uniprot_id"])
    if not rows:
        raise ReferenceError(f"{path} has no regions")
    return len(rows)


# commands ====
def parse_list(value: str | None, allowed: tuple[str, ...], name: str) -> tuple[str, ...]:
    if not value:
        return allowed
    chosen = tuple(v.strip() for v in value.split(",") if v.strip())
    unknown = [v for v in chosen if v not in allowed]
    if unknown:
        raise SystemExit(f"unknown {name}: {', '.join(unknown)} (choose from {', '.join(allowed)})")
    return chosen


def command_install(args: argparse.Namespace) -> int:
    root = reference_root(args.root)
    steps = parse_list(args.only, STEPS, "step")
    populations = parse_list(args.populations, POPULATIONS, "population")
    print(f"reference root: {root}\nsteps: {', '.join(steps)}" + (f"\npopulations: {', '.join(populations)}" if "panukb" in steps else ""))
    if args.dry_run:
        return 0
    root.mkdir(parents=True, exist_ok=True)
    if "genome" in steps:
        install_genome(root, args.force)
    if "annotation" in steps:
        install_annotation(root, args.force, Path(args.usefuldata_rda) if args.usefuldata_rda else None, args.rscript)
    if "panukb" in steps:
        install_panukb(root, populations, args.force)
    return command_check(args)


def command_check(args: argparse.Namespace) -> int:
    root = reference_root(args.root)
    steps = parse_list(args.only, STEPS, "step")
    populations = parse_list(args.populations, POPULATIONS, "population")
    problems = 0

    def report(label: str, action) -> None:
        nonlocal problems
        try:
            detail = action()
            print(f"ok       {label}" + (f"  {detail}" if detail else ""))
        except (ReferenceError, OSError, ValueError, KeyError) as error:
            problems += 1
            print(f"MISSING  {label}  {error}")

    print(f"reference root: {root}")
    if "genome" in steps:
        report("genome/hg38ToHg19.over.chain.gz", lambda: check_chain(root / "genome" / "hg38ToHg19.over.chain.gz"))
        report("genome/GRCh37.fa(.fai)", lambda: check_fasta(root / "genome" / "GRCh37.fa"))
    if "annotation" in steps:
        report("annotation/gene_regions_GRCh38.tsv", lambda: f"{check_annotation(root):,} gene regions")
    if "panukb" in steps:
        for population in populations:
            report(f"panukb/UKBB.{population}.variant_index.parquet", lambda p=population: f"{check_panukb(root, p):,} variants")
    manifest = read_manifest(root)
    for relative, item in sorted(manifest["items"].items()):
        path = root / relative
        if path.is_file() and path.stat().st_size != item.get("bytes"):
            problems += 1
            print(f"CHANGED  {relative}  size {path.stat().st_size} differs from the recorded {item.get('bytes')}")
    print("all reference data are in place" if problems == 0 else f"{problems} problem(s): run install again")
    return 0 if problems == 0 else 1


def command_paths(args: argparse.Namespace) -> int:
    root = reference_root(args.root)
    paths = {"chain": root / "genome" / "hg38ToHg19.over.chain.gz", "grch37_fasta": root / "genome" / "GRCh37.fa",
             "grch37_fai": root / "genome" / "GRCh37.fa.fai", "gene_annotation": root / "annotation" / "gene_regions_GRCh38.tsv",
             "panukb_dir": root / "panukb"}
    if args.nextflow:
        print(" ".join(f"--{key} {value}" for key, value in paths.items()))
    else:
        for key, value in paths.items():
            print(f"{key}\t{value}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("install", "check", "paths"):
        p = sub.add_parser(name)
        p.add_argument("--root", help="reference root (default GWASFM_REFERENCE_ROOT or ~/.cache/gwasfinemapping/reference)")
        if name != "paths":
            p.add_argument("--only", help=f"comma-separated steps of {', '.join(STEPS)} (default all)")
            p.add_argument("--populations", help=f"comma-separated Pan-UKB populations (default {', '.join(POPULATIONS)})")
        if name == "install":
            p.add_argument("--force", action="store_true", help="reinstall what is already there")
            p.add_argument("--dry-run", action="store_true", help="show the root, steps and populations only")
            p.add_argument("--usefuldata-rda", help=f"a local {USEFULDATA_DATASET}.rda instead of downloading it")
            p.add_argument("--rscript", help="Rscript to read the .rda with")
        if name == "paths":
            p.add_argument("--nextflow", action="store_true", help="print as --chain ... --grch37_fasta ... options for nextflow run")
    args = parser.parse_args(argv)
    try:
        return {"install": command_install, "check": command_check, "paths": command_paths}[args.command](args)
    except (ReferenceError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
