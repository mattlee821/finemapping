#!/usr/bin/env python3
"""Discover GWASprocessing outputs and run a lightweight indexed-VCF preflight."""

import argparse
import csv
import gzip
import io
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Any


MANIFEST_COLUMNS = ("study_id", "vcf", "index", "metadata", "qc")
SUMMARY_COLUMNS = (
    "study_id",
    "ancestry",
    "genome_build",
    "sample_id",
    "index_record_count",
    "vcf_bytes",
    "index_bytes",
    "phenotype",
)


def discover_inputs(input_dir: str, max_inputs: int = 0) -> tuple[list[dict[str, str]], dict[str, Any]]:
    root = Path(input_dir).expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"Input directory does not exist: {root}")
    if max_inputs < 0:
        raise ValueError("max_inputs must be zero (unlimited) or a positive integer.")

    directories = sorted(
        (path for path in root.iterdir() if path.is_dir() and not path.name.startswith(".")),
        key=lambda path: path.name,
    )
    if not directories:
        raise ValueError(f"No study directories found under: {root}")

    selected = directories if max_inputs == 0 else directories[:max_inputs]
    entries = []
    errors = []
    for directory in selected:
        study_id = directory.name
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", study_id):
            errors.append(f"{study_id}: unsupported directory name")
            continue
        paths = {
            "study_id": study_id,
            "vcf": str(directory / f"{study_id}.gwas.vcf.gz"),
            "index": str(directory / f"{study_id}.gwas.vcf.gz.tbi"),
            "metadata": str(directory / f"{study_id}.metadata.json"),
            "qc": str(directory / f"{study_id}.qc.json"),
        }
        missing = [key for key in ("vcf", "index", "metadata", "qc") if not Path(paths[key]).is_file()]
        if missing:
            errors.append(f"{study_id}: missing " + ", ".join(missing))
            continue
        entries.append(paths)

    if errors:
        raise ValueError("Invalid selected study directories:\n  " + "\n  ".join(errors))
    if not entries:
        raise ValueError(f"No valid study directories selected under: {root}")

    summary = {
        "input_dir": str(root),
        "discovered_count": len(directories),
        "selected_count": len(entries),
        "max_inputs": max_inputs,
        "selected_studies": [entry["study_id"] for entry in entries],
    }
    return entries, summary


def write_manifest(
    input_dir: str,
    output_path: str,
    summary_path: str,
    max_inputs: int = 0,
    existing_manifest: str | None = None,
    existing_summary: str | None = None,
) -> dict[str, Any]:
    entries, summary = discover_inputs(input_dir, max_inputs)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=MANIFEST_COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(entries)
    manifest_text = buffer.getvalue()
    if existing_manifest and Path(existing_manifest).is_file():
        if Path(existing_manifest).read_bytes() != manifest_text.encode("utf-8"):
            raise ValueError(
                "The existing output manifest differs from the selected inputs; "
                "use a new output directory."
            )
    if existing_summary and Path(existing_summary).is_file():
        existing_value = json.loads(Path(existing_summary).read_text(encoding="utf-8"))
        if existing_value != summary:
            raise ValueError(
                "The existing discovery summary differs from the selected inputs; "
                "use a new output directory."
            )
    output.write_text(manifest_text, encoding="utf-8", newline="")
    write_json(summary_path, summary)
    print(
        f"Discovered {summary['discovered_count']} study directories; "
        f"selected {summary['selected_count']} for this run."
    )
    return summary


def normalized_build(value: Any) -> str:
    normalized = re.sub(r"[^a-z0-9]", "", str(value).lower())
    if normalized.startswith(("grch38", "hg38")):
        return "GRCh38"
    if normalized.startswith(("grch37", "hg19")):
        return "GRCh37"
    return ""


def read_vcf_header(vcf_path: Path) -> tuple[list[str], str]:
    reference_values = []
    with gzip.open(vcf_path, "rt", encoding="utf-8") as stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if line.startswith("##reference="):
                reference_values.append(line.split("=", 1)[1])
            elif line.startswith("#CHROM\t"):
                columns = line[1:].split("\t")
                if len(reference_values) != 1:
                    raise ValueError(f"Expected exactly one ##reference declaration in {vcf_path}")
                fixed = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                if columns[: len(fixed)] != fixed or len(columns) != len(fixed) + 1:
                    raise ValueError(f"Expected standard GWAS-VCF columns and exactly one sample in {vcf_path}")
                return columns, reference_values[0]
            elif not line.startswith("#"):
                raise ValueError(f"VCF data appeared before the #CHROM header in {vcf_path}")
    raise ValueError(f"VCF header is missing the #CHROM column line: {vcf_path}")


def indexed_record_count(vcf_path: Path) -> int:
    try:
        result = subprocess.run(
            ["bcftools", "index", "--nrecords", str(vcf_path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as error:
        raise ValueError("bcftools is required for indexed-VCF preflight.") from error
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip() or error.stdout.strip()
        raise ValueError(f"bcftools could not read the VCF index for {vcf_path}: {detail}") from error
    try:
        count = int(result.stdout.strip())
    except ValueError as error:
        raise ValueError(f"bcftools returned an invalid record count for {vcf_path}") from error
    if count < 1:
        raise ValueError(f"VCF index contains no records: {vcf_path}")
    return count


def preflight_study(
    study_id: str,
    vcf_path: str,
    index_path: str,
    metadata_path: str,
    qc_path: str,
    output_path: str,
    source_vcf_path: str | None = None,
    source_index_path: str | None = None,
    source_metadata_path: str | None = None,
    source_qc_path: str | None = None,
) -> dict[str, Any]:
    vcf = Path(vcf_path)
    index = Path(index_path)
    if not vcf.is_file() or not index.is_file():
        raise ValueError(f"VCF or tabix index is missing for study {study_id}")

    metadata = json.loads(Path(metadata_path).read_text(encoding="utf-8"))
    qc = json.loads(Path(qc_path).read_text(encoding="utf-8"))
    if not isinstance(metadata, dict) or not isinstance(qc, dict):
        raise ValueError(f"Metadata and QC files must contain JSON objects for study {study_id}")
    if metadata.get("output_id") != study_id:
        raise ValueError(f"Metadata output_id does not match study directory: {study_id}")
    if metadata.get("standard_format") != "GWAS-VCF":
        raise ValueError(f"Metadata does not identify a GWAS-VCF for study {study_id}")
    ancestry = metadata.get("ancestry")
    if not isinstance(ancestry, str) or not ancestry.strip():
        raise ValueError(f"Metadata ancestry is missing for study {study_id}")

    build_detection = qc.get("build_detection")
    if not isinstance(build_detection, dict):
        raise ValueError(f"QC build_detection is missing for study {study_id}")
    builds = (
        metadata.get("coordinate_build"),
        qc.get("coordinate_build"),
        build_detection.get("build"),
    )
    if any(normalized_build(build) != "GRCh38" for build in builds):
        raise ValueError(f"Metadata/QC do not consistently identify GRCh38 for study {study_id}")

    columns, header_reference = read_vcf_header(vcf)
    if normalized_build(header_reference) != "GRCh38":
        raise ValueError(f"VCF reference does not identify GRCh38 for study {study_id}")

    count = indexed_record_count(vcf)
    report = {
        "schema_version": 1,
        "study_id": study_id,
        "ancestry": ancestry,
        "genome_build": "GRCh38",
        "vcf": str(Path(source_vcf_path or vcf_path).resolve()),
        "vcf_bytes": vcf.stat().st_size,
        "index": str(Path(source_index_path or index_path).resolve()),
        "index_bytes": index.stat().st_size,
        "index_record_count": count,
        "sample_id": columns[-1],
        "header_reference": header_reference,
        "phenotype": metadata.get("phenotype"),
        "metadata_file": str(Path(source_metadata_path or metadata_path).resolve()),
        "qc_file": str(Path(source_qc_path or qc_path).resolve()),
    }
    write_json(output_path, report)
    return report


def write_summary(report_paths: list[str], output_path: str) -> None:
    reports = [json.loads(Path(path).read_text(encoding="utf-8")) for path in report_paths]
    if not reports:
        raise ValueError("No preflight reports were provided.")
    study_ids = [report["study_id"] for report in reports]
    if len(study_ids) != len(set(study_ids)):
        raise ValueError("Duplicate study IDs found in preflight reports.")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=SUMMARY_COLUMNS, delimiter="\t")
        writer.writeheader()
        for report in sorted(reports, key=lambda value: value["study_id"]):
            writer.writerow({column: report.get(column, "") for column in SUMMARY_COLUMNS})


def write_json(path: str, value: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_path = Path(stream.name)
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary_path, output)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    manifest_parser = subparsers.add_parser("manifest")
    manifest_parser.add_argument("--input-dir", required=True)
    manifest_parser.add_argument("--output", required=True)
    manifest_parser.add_argument("--summary", required=True)
    manifest_parser.add_argument("--existing-manifest")
    manifest_parser.add_argument("--existing-summary")
    manifest_parser.add_argument("--max-inputs", type=int, default=0)

    preflight_parser = subparsers.add_parser("preflight")
    preflight_parser.add_argument("--study-id", required=True)
    preflight_parser.add_argument("--vcf", required=True)
    preflight_parser.add_argument("--index", required=True)
    preflight_parser.add_argument("--metadata", required=True)
    preflight_parser.add_argument("--qc", required=True)
    preflight_parser.add_argument("--source-vcf")
    preflight_parser.add_argument("--source-index")
    preflight_parser.add_argument("--source-metadata")
    preflight_parser.add_argument("--source-qc")
    preflight_parser.add_argument("--output", required=True)

    summary_parser = subparsers.add_parser("summarize")
    summary_parser.add_argument("--reports", nargs="+", required=True)
    summary_parser.add_argument("--output", required=True)

    args = parser.parse_args()
    try:
        if args.command == "manifest":
            write_manifest(
                args.input_dir,
                args.output,
                args.summary,
                args.max_inputs,
                args.existing_manifest,
                args.existing_summary,
            )
        elif args.command == "preflight":
            preflight_study(
                args.study_id,
                args.vcf,
                args.index,
                args.metadata,
                args.qc,
                args.output,
                args.source_vcf,
                args.source_index,
                args.source_metadata,
                args.source_qc,
            )
        else:
            write_summary(args.reports, args.output)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
