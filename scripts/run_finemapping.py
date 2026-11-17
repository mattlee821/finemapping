#!/usr/bin/env python3
"""Run a bounded, allele-aware Pan-UKB fine-mapping analysis for one GWAS-VCF."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


QUERY_FORMAT = (
    "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t"
    "[%LP{0}]\\t[%ES{0}]\\t[%SE{0}]\\t[%SS{0}]\\n"
)
METHODS = "p_value,p_ld,susie_rss,susie_rrss"
REFERENCE_SAMPLES = {"EUR": 420_542}


def resolve_reference_sample_size(population: str, override: int | None) -> int | None:
    if override is not None:
        if override < 3:
            raise ValueError("LD reference sample size must be at least 3.")
        return override
    return REFERENCE_SAMPLES.get(population)


def normalize_chromosome(chromosome: str) -> str:
    value = str(chromosome)
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def chromosome_sort_key(chromosome: str) -> tuple[int, str]:
    value = normalize_chromosome(chromosome)
    if value.isdigit():
        return int(value), ""
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(value.upper(), 26), value


def numeric(value: str) -> float | None:
    if value in ("", "."):
        return None
    result = float(value)
    if not math.isfinite(result):
        return None
    return result


def query_vcf(vcf: Path, region: str | None = None,
              minimum_lp: float | None = None) -> list[dict[str, Any]]:
    command = ["bcftools", "query"]
    if region:
        command.extend(["-r", region])
    if minimum_lp is not None:
        command.extend(["-i", f"FMT/LP[0:0]>={minimum_lp:.12g}"])
    command.extend(["-f", QUERY_FORMAT, str(vcf)])
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    records = []
    for line_number, line in enumerate(result.stdout.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(f"Unexpected bcftools query row {line_number}: {line}")
        chromosome, position, rsid, ref, alt, lp, effect, se, sample_size = fields
        records.append({
            "CHROM": normalize_chromosome(chromosome),
            "POS": int(position),
            "ID": rsid,
            "REF": ref,
            "ALT": alt,
            "LP": numeric(lp),
            "ES": numeric(effect),
            "SE": numeric(se),
            "SS": numeric(sample_size),
        })
    return records


def cluster_significant(records: list[dict[str, Any]], merge_distance_bp: int
                        ) -> list[list[dict[str, Any]]]:
    significant = [record for record in records if record["LP"] is not None]
    significant.sort(key=lambda r: (chromosome_sort_key(r["CHROM"]), r["POS"]))
    clusters: list[list[dict[str, Any]]] = []
    for record in significant:
        if (
            not clusters
            or clusters[-1][-1]["CHROM"] != record["CHROM"]
            or record["POS"] - clusters[-1][-1]["POS"] > merge_distance_bp
        ):
            clusters.append([record])
        else:
            clusters[-1].append(record)
    return clusters


def locus_window(cluster_start: int, cluster_end: int,
                 flank_bp: int) -> tuple[int, int]:
    if cluster_start < 1 or cluster_end < cluster_start or flank_bp < 1:
        raise ValueError("Cluster coordinates and flank must be positive and ordered.")
    return max(1, cluster_start - flank_bp), cluster_end + flank_bp


def select_mapped_loci(loci: list[dict[str, Any]],
                       significant: list[dict[str, Any]],
                       max_loci: int) -> list[dict[str, Any]]:
    if max_loci < 0:
        raise ValueError("Maximum locus count must be non-negative.")
    selected = []
    for locus in loci:
        if locus["liftover_status"] != "mapped":
            continue
        lead = next(
            row for row in significant
            if row["CHROM"] == locus["CHROM"]
            and row["POS"] == locus["lead_POS"]
            and row["ID"] == locus["lead_ID"]
        )
        selected.append({
            **lead,
            "cluster_start": locus["start"],
            "cluster_end": locus["end"],
        })
        if max_loci and len(selected) >= max_loci:
            break
    return selected


def map_variant(lifter: Any, record: dict[str, Any]
                ) -> tuple[dict[str, Any] | None, str | None]:
    ref = record["REF"].upper()
    alt = record["ALT"].upper()
    if not ref or not alt or re.search(r"[^ACGTN]", ref + alt):
        return None, "unsupported_symbolic_or_non_dna_allele"

    chromosome = normalize_chromosome(record["CHROM"])
    source_positions = range(record["POS"] - 1, record["POS"] - 1 + len(ref))
    mappings = [lifter.convert_coordinate(f"chr{chromosome}", p)
                for p in source_positions]
    if any(len(options) == 0 for options in mappings):
        return None, "no_chain_mapping"
    if any(len(options) != 1 for options in mappings):
        return None, "ambiguous_chain_mapping"

    mapped = [options[0] for options in mappings]
    target_chroms = {normalize_chromosome(item[0]) for item in mapped}
    strands = {item[2] for item in mapped}
    chain_ids = {item[3] for item in mapped}
    positions = [int(item[1]) for item in mapped]
    if len(target_chroms) != 1 or len(strands) != 1 or len(chain_ids) != 1:
        return None, "reference_span_crosses_chain_blocks"
    strand = next(iter(strands))
    step = 1 if strand == "+" else -1
    if any(b - a != step for a, b in zip(positions, positions[1:])):
        return None, "reference_span_not_contiguous_after_liftover"

    if strand == "-":
        complement = str.maketrans("ACGTN", "TGCAN")
        ref = ref.translate(complement)[::-1]
        alt = alt.translate(complement)[::-1]
    lifted = dict(record)
    lifted.update({
        "LD_CHROM": next(iter(target_chroms)),
        "LD_POS": min(positions) + 1,
        "LD_REF": ref,
        "LD_ALT": alt,
        "LIFTOVER_STRAND": strand,
        "LIFTOVER_CHAIN": next(iter(chain_ids)),
    })
    return lifted, None


def read_gene_annotation(path: Path, uniprot_id: str) -> dict[str, Any]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    matches = [row for row in rows if row["uniprot_id"] == uniprot_id]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one gene annotation for {uniprot_id}; found {len(matches)}"
        )
    gene = matches[0]
    if gene["genome_build"] != "GRCh38":
        raise ValueError("Gene annotation must use GRCh38 coordinates.")
    gene["chrom"] = normalize_chromosome(gene["chrom"])
    gene["start"] = int(gene["start"])
    gene["end"] = int(gene["end"])
    if gene["start"] < 1 or gene["end"] < gene["start"]:
        raise ValueError(f"Invalid gene coordinates for {uniprot_id}")
    return gene


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if row.get(key) is None else row.get(key)
                             for key in columns})


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_checked(command: list[str], env: dict[str, str] | None = None,
                cwd: Path | None = None
                ) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        command, capture_output=True, text=True, env=env, cwd=cwd
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: "
            f"{' '.join(command)}\n{result.stderr.strip()}"
        )
    return result


def write_lifted_vcf(records: list[dict[str, Any]], lifter: Any,
                     output_path: Path) -> tuple[dict[str, dict[str, Any]],
                                                list[dict[str, Any]]]:
    lifted_by_id: dict[str, dict[str, Any]] = {}
    excluded = []
    eligible = []
    for i, record in enumerate(records):
        unique_id = f"v{i:08d}"
        lifted, reason = map_variant(lifter, record)
        if reason:
            excluded.append({**record, "status": "excluded", "reason": reason})
            continue
        values = [record[field] for field in ("ES", "SE", "LP", "SS")]
        if any(value is None for value in values[:3]):
            excluded.append({
                **record,
                "status": "excluded",
                "reason": "missing_ES_SE_or_LP",
            })
            continue
        lifted_by_id[unique_id] = lifted
        eligible.append((unique_id, record, lifted, values))

    with output_path.open("w", encoding="utf-8") as stream:
        stream.write(
            "##fileformat=VCFv4.2\n"
            "##FORMAT=<ID=ES,Number=1,Type=Float,Description=\"Effect size\">\n"
            "##FORMAT=<ID=SE,Number=1,Type=Float,Description=\"Standard error\">\n"
            "##FORMAT=<ID=LP,Number=1,Type=Float,Description=\"-log10 p-value\">\n"
            "##FORMAT=<ID=SS,Number=1,Type=Float,Description=\"Sample size\">\n"
        )
        for chromosome in sorted({
            lifted["LD_CHROM"] for _, _, lifted, _ in eligible
        }, key=chromosome_sort_key):
            stream.write(f"##contig=<ID={chromosome}>\n")
        stream.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGWAS\n"
        )
        for unique_id, record, lifted, values in eligible:
            sample = ":".join(str(value) if value is not None else "." for value in values)
            stream.write(
                f"{lifted['LD_CHROM']}\t{lifted['LD_POS']}\t{unique_id}\t"
                f"{lifted['LD_REF']}\t{lifted['LD_ALT']}\t.\tPASS\t.\t"
                f"ES:SE:LP:SS\t{sample}\n"
            )
    return lifted_by_id, excluded


def normalize_lifted_vcf(vcf: Path, fasta: Path, output_path: Path
                         ) -> list[dict[str, Any]]:
    with tempfile.TemporaryDirectory(prefix="gwasfm-normalize-") as directory:
        normalized = Path(directory) / "normalized.vcf"
        run_checked([
            "bcftools", "norm", "-f", str(fasta), "-c", "w", "-Ov",
            "-o", str(normalized), str(vcf),
        ])
        result = run_checked([
            "bcftools", "query", "-f",
            "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t"
            "[%ES{0}]\\t[%SE{0}]\\t[%LP{0}]\\t[%SS{0}]\\n",
            str(normalized),
        ])
    rows = []
    for line in result.stdout.splitlines():
        chrom, pos, unique_id, ref, alt, es, se, lp, ss = line.split("\t")
        rows.append({
            "source_variant_key": unique_id,
            "LD_CHROM": normalize_chromosome(chrom),
            "LD_POS": int(pos),
            "LD_REF": ref,
            "LD_ALT": alt,
            "ES": numeric(es),
            "SE": numeric(se),
            "LP": numeric(lp),
            "SS": numeric(ss),
        })
    write_tsv(
        output_path,
        rows,
        ["source_variant_key", "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT",
         "ES", "SE", "LP", "SS"],
    )
    return rows


def resolve_panel_matches(normalized: list[dict[str, Any]], matches: Any
                          ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    grouped: dict[int, list[dict[str, Any]]] = {}
    excluded = []
    for row, match in zip(normalized, matches.itertuples(index=False)):
        panel_index = int(match.idx)
        if panel_index < 0:
            excluded.append({
                **row, "status": "excluded", "reason": "no_exact_panel_allele_match"
            })
            continue
        candidate = {
            **row,
            "PANEL_IDX": panel_index,
            "PANEL_ALLELE_SWAPPED": bool(match.flip),
            "LD_VARIANT_ID": (
                f"{row['LD_CHROM']}:{row['LD_POS']}:{row['LD_REF']}:{row['LD_ALT']}"
            ),
        }
        grouped.setdefault(panel_index, []).append(candidate)

    matched = []
    for panel_index, candidates in grouped.items():
        candidates.sort(
            key=lambda row: (
                -(row["LP"] if row["LP"] is not None else float("-inf")),
                row["source_variant_key"],
            )
        )
        winner = candidates[0]
        matched.append(winner)
        for duplicate in candidates[1:]:
            excluded.append({
                **duplicate,
                "status": "excluded",
                "reason": "duplicate_panel_variant",
                "duplicate_of": winner["source_variant_key"],
            })
    return matched, excluded


def default_panukb_dir() -> Path:
    """The Pan-UKB folder of the reference root that scripts/install_reference_data.sh fills."""
    root = os.environ.get("GWASFM_REFERENCE_ROOT") or "~/.cache/gwasfinemapping/reference"
    return Path(root).expanduser() / "panukb"


def installed_population_index(population: str, panukb_dir: Path | None) -> Path | None:
    """The whole-population variant index of the reference root, if it is installed and its record says it is this
    population on GRCh37 for the whole genome; None otherwise (the regional build below is then used)."""
    if panukb_dir is None:
        return None
    index = panukb_dir / f"UKBB.{population}.variant_index.parquet"
    sidecar = Path(f"{index}.json")
    if not index.is_file() or not sidecar.is_file():
        return None
    record = json.loads(sidecar.read_text())
    if (record.get("population"), record.get("region"), record.get("build")) != (population, "genome", "GRCh37"):
        raise ValueError(f"{sidecar} does not describe the whole-genome {population} index on GRCh37")
    return index


def cut_regional_index(whole: Path, chromosome: str, start: int, end: int, destination: Path) -> int:
    """Write the rows of a whole-population index inside chromosome:start-end to destination (same columns, same
    types), reading only the matching row groups. Returns the row count."""
    import pyarrow.parquet as pq

    table = pq.read_table(whole, filters=[("contig", "=", chromosome), ("position", ">=", start),
                                          ("position", "<=", end)])
    temporary = destination.with_name(destination.name + ".tmp")
    pq.write_table(table, temporary)
    os.replace(temporary, destination)
    return table.num_rows


def ensure_regional_index(population: str, chromosome: str, start: int, end: int,
                          cache_dir: Path, scripts_dir: Path,
                          panukb_dir: Path | None = None) -> Path:
    """The variant index of the LD panel for 1 region. It is cut from the installed whole-population index when there
    is one (no Hail, no network); otherwise it is built from the public Pan-UKB table through Hail and S3."""
    region = f"{chromosome}:{start}-{end}"
    index = cache_dir / "panukb" / (
        f"UKBB.{population}.chr{chromosome}_{start}_{end}.variant_index.parquet"
    )
    sidecar = Path(f"{index}.json")
    expected = {"population": population, "region": region, "build": "GRCh37"}
    index.parent.mkdir(parents=True, exist_ok=True)
    if index.exists():
        if not sidecar.is_file() or json.loads(sidecar.read_text()) != expected:
            raise ValueError(f"Existing variant index lacks matching provenance: {index}")
        return index
    whole = installed_population_index(population, panukb_dir)
    if whole is not None:
        rows = cut_regional_index(whole, chromosome, start, end, index)
        if rows == 0:
            index.unlink()
            raise ValueError(f"The installed {population} index has no variants in {region}")
    else:
        run_checked([
            "bash", str(scripts_dir / "build_panukb_variant_index.sh"),
            population, str(index), region,
        ])
    write_json(sidecar, expected)
    return index


def matrix_to_tsv(npz_path: Path, variants_path: Path, output_path: Path
                  ) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    ld_variants = pd.read_csv(variants_path, sep="\t", dtype={"chrom": str})
    with np.load(npz_path) as data:
        matrix = data["matrix"]
        indexes = data["idx"]
    if matrix.shape != (len(ld_variants), len(ld_variants)):
        raise ValueError("Extracted LD dimensions do not match its variant table.")
    if len(indexes) != len(ld_variants) or np.any(~np.isfinite(matrix)):
        raise ValueError("Extracted LD has missing values or inconsistent indices.")
    ids = [
        f"{normalize_chromosome(row.chrom)}:{int(row.pos)}:{row.ref}:{row.alt}"
        for row in ld_variants.itertuples(index=False)
    ]
    if len(set(ids)) != len(ids):
        raise ValueError("Extracted LD contains duplicate variant identifiers.")
    if np.max(np.abs(matrix - matrix.T)) > 1e-8:
        raise ValueError("Extracted LD matrix is not symmetric.")
    diagonal = np.diag(matrix)
    if np.any(diagonal <= 0) or np.max(np.abs(diagonal - 1)) > 0.01:
        raise ValueError("Pan-UKB LD diagonal is outside the 0.01 normalization tolerance.")
    scale = np.sqrt(diagonal)
    normalized = matrix / scale[:, None] / scale[None, :]
    np.fill_diagonal(normalized, 1.0)
    if np.any(~np.isfinite(normalized)) or np.max(np.abs(normalized)) > 1.00001:
        raise ValueError("Normalized LD contains invalid correlations.")

    normalized_npz = output_path.with_name("normalized_matrix.npz")
    np.savez_compressed(normalized_npz, matrix=normalized, idx=indexes)
    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["LD_VARIANT_ID", *ids])
        for variant_id, row in zip(ids, normalized):
            writer.writerow([variant_id, *(f"{value:.12g}" for value in row)])
    return {
        "variants": len(ids),
        "shape": list(matrix.shape),
        "raw_diagonal_min": float(diagonal.min()),
        "raw_diagonal_max": float(diagonal.max()),
        "diagonal_normalization": "D^-1/2 R D^-1/2",
        "max_abs_asymmetry": float(np.max(np.abs(matrix - matrix.T))),
        "missing_cells": 0,
    }


def run_methods(variants_path: Path, matrix_path: Path | None, output_dir: Path,
                reference_n: int | None, region_id: str, region_type: str,
                methods: str = METHODS) -> None:
    project_root = Path(__file__).resolve().parents[1]
    command = [
        "Rscript", str(Path(__file__).with_name("run_region_methods.R")),
        "--variants", str(variants_path.resolve()),
        "--outdir", str(output_dir.resolve()),
        "--methods", methods,
        "--region-id", region_id,
        "--region-type", region_type,
    ]
    if reference_n is not None:
        command.extend(["--reference-sample-size", str(reference_n)])
    if matrix_path is not None:
        command.extend(["--matrix", str(matrix_path.resolve())])
    environment = os.environ.copy()
    tools_root = Path(
        environment.get(
            "GWASFM_TOOLS_ROOT", str(Path.home() / ".cache" / "gwasfinemapping")
        )
    ).expanduser()
    if not tools_root.is_absolute():
        tools_root = project_root / tools_root
    environment["RENV_PATHS_LIBRARY"] = str(tools_root / "renv-library")
    run_checked(command, env=environment, cwd=project_root)


def aggregate_method_results(output_root: Path, scripts_dir: Path) -> None:
    project_root = scripts_dir.parent
    command = [
        "Rscript", str(scripts_dir / "aggregate_method_results.R"),
        "--output-root", str(output_root.resolve()),
    ]
    environment = os.environ.copy()
    tools_root = Path(
        environment.get(
            "GWASFM_TOOLS_ROOT", str(Path.home() / ".cache" / "gwasfinemapping")
        )
    ).expanduser()
    if not tools_root.is_absolute():
        tools_root = project_root / tools_root
    environment["RENV_PATHS_LIBRARY"] = str(tools_root / "renv-library")
    run_checked(command, env=environment, cwd=project_root)


def analyze_locus(locus: dict[str, Any], vcf: Path, lifter: Any, fasta: Path,
                  population: str, reference_sample_size: int, cache_dir: Path,
                  output_dir: Path, scripts_dir: Path, flank_bp: int,
                  gene: dict[str, Any], panukb_dir: Path | None = None) -> dict[str, Any]:
    chrom = locus["CHROM"]
    start, end = locus_window(
        locus["cluster_start"], locus["cluster_end"], flank_bp
    )
    region = f"{chrom}:{start}-{end}"
    records = query_vcf(vcf, region=region)
    if not records:
        raise ValueError(f"No GWAS variants found in selected locus {region}")

    locus_dir = output_dir / "loci" / (
        f"chr{chrom}_{start}_{end}"
    )
    inputs_dir = locus_dir / "inputs"
    mapping_dir = locus_dir / "mapping"
    ld_dir = locus_dir / "ld"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    mapping_dir.mkdir(parents=True, exist_ok=True)
    ld_dir.mkdir(parents=True, exist_ok=True)

    lifted_vcf = inputs_dir / "lifted.vcf"
    lifted_by_id, excluded = write_lifted_vcf(records, lifter, lifted_vcf)
    normalized_path = mapping_dir / "normalized_grch37_variants.tsv"
    normalized = normalize_lifted_vcf(lifted_vcf, fasta, normalized_path)
    normalized = [
        {**lifted_by_id[row["source_variant_key"]], **row}
        for row in normalized
        if row["source_variant_key"] in lifted_by_id
    ]
    if not normalized:
        raise ValueError(f"No variants mapped to GRCh37 in selected locus {region}")

    target_chromosomes = {row["LD_CHROM"] for row in normalized}
    if len(target_chromosomes) != 1:
        raise ValueError("A fine-mapping locus mapped to multiple chromosomes.")
    target_chrom = next(iter(target_chromosomes))
    target_start = max(1, min(row["LD_POS"] for row in normalized) - 100)
    target_end = max(row["LD_POS"] for row in normalized) + 100
    index = ensure_regional_index(
        population, target_chrom, target_start, target_end, cache_dir, scripts_dir, panukb_dir
    )

    import pandas as pd
    from ldcov.io.variant_index import VariantIndex

    candidates = pd.DataFrame({
        "contig": [row["LD_CHROM"] for row in normalized],
        "position": [row["LD_POS"] for row in normalized],
        "ref": [row["LD_REF"] for row in normalized],
        "alt": [row["LD_ALT"] for row in normalized],
    })
    matches = VariantIndex(str(index)).match_variants(
        candidates, allow_allele_swap=True
    )
    matched, panel_excluded = resolve_panel_matches(normalized, matches)
    if not matched:
        raise ValueError(f"No variants matched the {population} LD panel in {region}")

    matched.sort(key=lambda r: (chromosome_sort_key(r["LD_CHROM"]), r["LD_POS"],
                                r["LD_REF"], r["LD_ALT"]))
    original_by_key = {
        (row["CHROM"], row["POS"], row["REF"], row["ALT"]): row for row in records
    }
    matched_rows = []
    z_path = inputs_dir / "panukb.z"
    with z_path.open("w", encoding="utf-8") as stream:
        stream.write("rsid chromosome position allele1 allele2\n")
        for i, row in enumerate(matched):
            source = original_by_key.get(
                (row["CHROM"], row["POS"], row["REF"], row["ALT"])
            )
            if source is None:
                raise ValueError("Could not recover the source GWAS variant row.")
            out = {
                **source,
                "GWAS_CHROM": source["CHROM"],
                "GWAS_POS": source["POS"],
                "GWAS_REF": source["REF"],
                "GWAS_ALT": source["ALT"],
                "CHROM": row["LD_CHROM"],
                "POS": row["LD_POS"],
                "REF": row["LD_REF"],
                "ALT": row["LD_ALT"],
                "LD_CHROM": row["LD_CHROM"],
                "LD_POS": row["LD_POS"],
                "LD_REF": row["LD_REF"],
                "LD_ALT": row["LD_ALT"],
                "LD_VARIANT_ID": row["LD_VARIANT_ID"],
                "PANEL_IDX": row["PANEL_IDX"],
                "PANEL_ALLELE_SWAPPED": row["PANEL_ALLELE_SWAPPED"],
                "EZ": source["ES"] / source["SE"],
                "cis_trans": (
                    "cis"
                    if source["CHROM"] == gene["chrom"]
                    and gene["start"] - 1_000_000 <= source["POS"]
                    <= gene["end"] + 1_000_000
                    else "trans"
                ),
            }
            matched_rows.append(out)
            stream.write(
                f"v{i:08d} {row['LD_CHROM']} {row['LD_POS']} "
                f"{row['LD_REF']} {row['LD_ALT']}\n"
            )

    mapping_rows = [
        {**row, "status": "matched", "reason": ""}
        for row in matched
    ] + excluded + panel_excluded
    mapping_columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "LP", "ES", "SE", "SS",
        "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "LIFTOVER_STRAND",
        "LIFTOVER_CHAIN", "PANEL_IDX", "PANEL_ALLELE_SWAPPED", "status", "reason",
        "duplicate_of",
    ]
    write_tsv(mapping_dir / "variant_mapping.tsv", mapping_rows, mapping_columns)
    write_tsv(
        inputs_dir / "matched_variants.tsv",
        matched_rows,
        ["CHROM", "POS", "ID", "REF", "ALT", "GWAS_CHROM", "GWAS_POS",
         "GWAS_REF", "GWAS_ALT", "LP", "ES", "SE", "SS",
         "EZ", "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "LD_VARIANT_ID",
         "PANEL_IDX", "PANEL_ALLELE_SWAPPED", "cis_trans"],
    )

    ld_prefix = ld_dir / f"panukb_{population.lower()}"
    run_checked([
        "ldcov", "--ld-bm", "--bm",
        f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.bm",
        "--variant-index", str(index),
        "--z", str(z_path),
        "--allow-allele-swap",
        "--on-missing", "error",
        "--block-cache", "2",
        "--fetch-workers", "2",
        "--output-format", "npz",
        "--out", str(ld_prefix),
    ])
    matrix_path = ld_dir / "matrix.tsv.gz"
    validation = matrix_to_tsv(
        Path(f"{ld_prefix}.npz"),
        Path(f"{ld_prefix}.variants.tsv"),
        matrix_path,
    )
    write_json(ld_dir / "validation.json", validation)
    methods_dir = locus_dir / "methods"
    reference_n = reference_sample_size
    if reference_n is None:
        raise ValueError(
            f"Supply --reference-sample-size for the {args.population} LD panel."
        )
    run_methods(
        inputs_dir / "matched_variants.tsv",
        matrix_path,
        methods_dir,
        reference_n,
        region_id=locus_dir.name,
        region_type="genome_wide_locus",
    )
    return {
        "region_id": locus_dir.name,
        "source_region": region,
        "variants_in_source_region": len(records),
        "variants_lifted": len(normalized),
        "variants_ld_matched": len(matched_rows),
        "variants_excluded": len(records) - len(matched_rows),
        "lead_variant_id": locus["ID"],
        "lead_lp": locus["LP"],
        "cluster_start": locus["cluster_start"],
        "cluster_end": locus["cluster_end"],
        "ld_index": str(index),
    }


def run(args: argparse.Namespace) -> Path:
    vcf = Path(args.vcf).resolve()
    metadata_path = Path(args.metadata).resolve()
    qc_path = Path(args.qc).resolve()
    chain = Path(args.chain).resolve()
    fasta = Path(args.grch37_fasta).resolve()
    fasta_index = Path(f"{fasta}.fai")
    gene_path = Path(args.gene_annotation).resolve()
    scripts_dir = Path(__file__).resolve().parent
    if not all(path.is_file() for path in (vcf, metadata_path, qc_path, chain,
                                           fasta, fasta_index, gene_path)):
        raise FileNotFoundError("VCF, metadata, QC, liftover, and FASTA inputs must exist.")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    qc = json.loads(qc_path.read_text(encoding="utf-8"))
    gwas_id = metadata["output_id"]
    study_id = metadata["study_id"]
    if metadata.get("coordinate_build") != "GRCh38" or qc.get("coordinate_build") != "GRCh38":
        raise ValueError("The real-data runner requires a consistently declared GRCh38 GWAS-VCF.")
    uniprot_id = re.match(r"[A-Z0-9]+", gwas_id)
    if not uniprot_id:
        raise ValueError(f"Could not identify a UniProt accession from {gwas_id}")
    gene = read_gene_annotation(gene_path, uniprot_id.group(0))
    reference_n = resolve_reference_sample_size(args.population, args.reference_sample_size)
    if "susie_rrss" in METHODS.split(",") and reference_n is None:
        raise ValueError(
            f"Supply --reference-sample-size for the {args.population} LD panel."
        )
    gwas_ancestry = metadata.get("ancestry")
    if gwas_ancestry in REFERENCE_SAMPLES and gwas_ancestry != args.population:
        raise ValueError(
            f"GWAS ancestry {gwas_ancestry} does not match LD population {args.population}."
        )
    if gwas_ancestry == "ALL":
        raise ValueError(
            "Pan-UKB has no combined-ancestry LD panel; use a matched reference or "
            "P-value-only analysis for ALL-ancestry GWAS."
        )
    output_dir = Path(args.output_root).resolve() / study_id / gwas_id
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"GWAS output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
    significant = query_vcf(vcf, minimum_lp=args.significance_lp)
    significant.sort(key=lambda r: (-r["LP"], chromosome_sort_key(r["CHROM"]), r["POS"]))
    gene_start = max(1, gene["start"] - args.gene_window_bp)
    gene_end = gene["end"] + args.gene_window_bp
    for row in significant:
        row["cis_trans"] = (
            "cis" if row["CHROM"] == gene["chrom"] and gene_start <= row["POS"] <= gene_end
            else "trans"
        )
    significant_columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "LP", "ES", "SE", "SS", "cis_trans"
    ]
    write_tsv(output_dir / "discovery" / "genomewide_significant.tsv",
              significant, significant_columns)

    gene_region = query_vcf(
        vcf, region=f"{gene['chrom']}:{gene_start}-{gene_end}"
    )
    for row in gene_region:
        row["EZ"] = row["ES"] / row["SE"] if row["ES"] is not None and row["SE"] else None
        row["cis_trans"] = "cis"
    gene_variants_path = output_dir / "gene_region" / "variants.tsv"
    write_tsv(
        gene_variants_path,
        gene_region,
        ["CHROM", "POS", "ID", "REF", "ALT", "LP", "ES", "SE", "SS",
         "EZ", "cis_trans"],
    )
    run_methods(
        gene_variants_path,
        None,
        output_dir / "gene_region" / "methods",
        reference_n,
        region_id="gene_region",
        region_type="gene_region",
        methods="p_value",
    )

    clusters = cluster_significant(significant, args.merge_distance_bp)
    loci = []
    for cluster in clusters:
        lead = max(cluster, key=lambda r: (r["LP"], -r["POS"]))
        mapped, reason = map_variant(lifter, lead)
        loci.append({
            "CHROM": lead["CHROM"],
            "start": min(row["POS"] for row in cluster),
            "end": max(row["POS"] for row in cluster),
            "significant_variant_count": len(cluster),
            "lead_ID": lead["ID"],
            "lead_POS": lead["POS"],
            "lead_LP": lead["LP"],
            "lead_cis_trans": lead["cis_trans"],
            "lead_GRCh37_CHROM": mapped["LD_CHROM"] if mapped else "",
            "lead_GRCh37_POS": mapped["LD_POS"] if mapped else "",
            "liftover_status": "mapped" if mapped else "excluded",
            "liftover_reason": reason or "",
        })
    loci.sort(key=lambda row: (-row["lead_LP"],
                               chromosome_sort_key(row["CHROM"]), row["lead_POS"]))
    write_tsv(
        output_dir / "discovery" / "genomewide_loci.tsv",
        loci,
        ["CHROM", "start", "end", "significant_variant_count", "lead_ID",
         "lead_POS", "lead_LP", "lead_cis_trans", "lead_GRCh37_CHROM",
         "lead_GRCh37_POS", "liftover_status", "liftover_reason"],
    )

    selected_loci = select_mapped_loci(loci, significant, args.max_loci)
    if not selected_loci:
        raise ValueError("No genome-wide significant locus maps to the GRCh37 reference.")
    analyzed = []
    for locus in selected_loci:
        analyzed.append(analyze_locus(
            locus=locus,
            vcf=vcf,
            lifter=lifter,
            fasta=fasta,
            population=args.population,
            reference_sample_size=reference_n,
            cache_dir=Path(args.cache_dir).expanduser().resolve(),
            panukb_dir=(Path(args.panukb_dir).expanduser().resolve() if args.panukb_dir else default_panukb_dir()),
            output_dir=output_dir,
            scripts_dir=scripts_dir,
            flank_bp=args.region_flank_bp,
            gene=gene,
        ))
    aggregate_method_results(output_dir, scripts_dir)

    sample_sizes = [row["SS"] for row in significant if row["SS"] is not None]
    provenance = {
        "schema_version": 1,
        "study_id": study_id,
        "gwas_id": gwas_id,
        "trait": metadata.get("phenotype"),
        "gwas_ancestry": metadata.get("ancestry"),
        "gwas_build": "GRCh38",
        "ld_reference": "Pan-UKB public AWS",
        "ld_population": args.population,
        "ld_reference_sample_size": reference_n,
        "reference_sample_size_source": (
            "EUR.ldadj.variant.ht global n_samples" if args.population == "EUR"
            and args.reference_sample_size is None else "explicit run override"
        ),
        "liftover_chain": str(chain),
        "liftover_chain_sha256": sha256_file(chain),
        "target_reference_fasta": str(fasta),
        "target_reference_fai": str(fasta_index),
        "gene_annotation": gene,
        "cis_window_bp": args.gene_window_bp,
        "significance_threshold_p": 10 ** (-args.significance_lp),
        "significance_threshold_lp": args.significance_lp,
        "locus_merge_distance_bp": args.merge_distance_bp,
        "fine_mapping_window": (
            f"cluster start/end +/- {args.region_flank_bp} bp"
        ),
        "fine_mapping_window_definition": "significant_cluster_bounds_plus_flank",
        "max_loci": args.max_loci,
        "genomewide_significant_variants": len(significant),
        "genomewide_loci": len(loci),
        "analyzed_loci": analyzed,
        "minimum_variant_ss_over_significant_variants": (
            min(sample_sizes) if sample_sizes else None
        ),
        "fine_mapping_methods": METHODS.split(","),
        "duplicate_panel_variant_policy": (
            "retain the highest-LP GWAS row per Pan-UKB variant index; "
            "record all duplicates in variant_mapping.tsv"
        ),
    }
    write_json(output_dir / "provenance.json", provenance)
    write_json(output_dir / "metadata.json", metadata)
    write_json(output_dir / "qc.json", qc)
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--chain", required=True, help="GRCh38-to-GRCh37 UCSC chain")
    parser.add_argument("--grch37-fasta", required=True)
    parser.add_argument("--gene-annotation", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--cache-dir", default="~/.cache/gwasfinemapping")
    parser.add_argument("--panukb-dir", help="folder with the installed whole-population indexes UKBB.<POP>.variant_index.parquet "
                        "(default: panukb/ under GWASFM_REFERENCE_ROOT, or ~/.cache/gwasfinemapping/reference); a regional "
                        "index is cut from it, and built through Hail only when it is not installed")
    parser.add_argument("--population", choices=("AFR", "AMR", "CSA", "EAS", "EUR", "MID"),
                        default="EUR")
    parser.add_argument("--reference-sample-size", type=int)
    parser.add_argument("--significance-lp", type=float, default=-math.log10(5e-8))
    parser.add_argument("--merge-distance-bp", type=int, default=250_000)
    parser.add_argument("--region-flank-bp", type=int, default=100_000)
    parser.add_argument("--gene-window-bp", type=int, default=1_000_000)
    parser.add_argument("--max-loci", type=int, default=0,
                        help="Maximum mapped loci to fine-map; 0 means all.")
    args = parser.parse_args()
    if args.max_loci < 0 or args.region_flank_bp < 1 or args.gene_window_bp < 0:
        parser.error("Locus count must be non-negative and window sizes valid.")
    if args.merge_distance_bp < 0 or args.significance_lp <= 0:
        parser.error("Merge distance must be non-negative and significance LP positive.")
    output = run(args)
    print(output)


if __name__ == "__main__":
    main()
