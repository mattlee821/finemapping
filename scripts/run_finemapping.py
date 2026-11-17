#!/usr/bin/env python3
"""Run scalable, allele-aware Pan-UKB statistical fine-mapping for one GWAS-VCF."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


BASE_QUERY_FORMAT = (
    "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t"
    "[%LP{0}]\\t[%ES{0}]\\t[%SE{0}]\\t[%SS{0}]"
)
METHODS = "p_value,p_ld,susie_rss,susie_rrss"

# Exact sample sizes verified from s3://pan-ukb-us-east-1/ld_release/UKBB.<POP>.ldadj.variant.ht/globals
REFERENCE_SAMPLES: dict[str, int] = {
    "AFR": 6_637,
    "AMR": 982,
    "CSA": 8_876,
    "EAS": 2_709,
    "EUR": 420_542,
    "MID": 1_599,
}

# Extended HLA single-locus bounds on GRCh38 (Sun et al. 2023, 8.50 Mb < 10.0 Mb Pan-UKB LD band)
MHC_CHROM = "6"
MHC_START_GRCH38 = 25_500_000
MHC_END_GRCH38 = 34_000_000

# Authoritative NCBI Genome Reference Consortium core MHC bounds (assembly_regions.txt)
MHC_CORE_GRC_GRCH38 = (28_510_120, 33_480_577)
MHC_CORE_GRC_GRCH37 = (28_477_797, 33_448_354)


def resolve_reference_sample_size(population: str, override: int | None) -> int | None:
    if override is not None:
        if override < 3:
            raise ValueError("LD reference sample size must be at least 3.")
        return override
    return REFERENCE_SAMPLES.get(population)


def normalize_chromosome(chromosome: str) -> str:
    value = str(chromosome).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return value


def chromosome_sort_key(chromosome: str) -> tuple[int, str]:
    value = normalize_chromosome(chromosome)
    if value.isdigit():
        return int(value), ""
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(value.upper(), 26), value


def numeric(value: str | None) -> float | None:
    if value is None or value in ("", ".", "NA"):
        return None
    result = float(value)
    if not math.isfinite(result):
        return None
    return result


def clean_str(value: str | None) -> str | None:
    if value is None or value in ("", "."):
        return None
    return value


def is_mhc_grch38(chrom: str, pos: int,
                  mhc_start: int = MHC_START_GRCH38,
                  mhc_end: int = MHC_END_GRCH38) -> bool:
    return normalize_chromosome(chrom) == MHC_CHROM and mhc_start <= pos <= mhc_end


def is_mhc_core_grc_grch38(chrom: str, pos: int) -> bool:
    return (
        normalize_chromosome(chrom) == MHC_CHROM
        and MHC_CORE_GRC_GRCH38[0] <= pos <= MHC_CORE_GRC_GRCH38[1]
    )


def inspect_vcf_tags(vcf: Path) -> dict[str, bool]:
    result = subprocess.run(
        ["bcftools", "view", "-h", str(vcf)],
        check=True,
        capture_output=True,
        text=True,
    )
    header = result.stdout
    return {
        "FORMAT_AF": bool(re.search(r"^##FORMAT=<ID=AF[>,]", header, re.M)),
        "FORMAT_AFS": bool(re.search(r"^##FORMAT=<ID=AFS[>,]", header, re.M)),
        "INFO_HM": bool(re.search(r"^##INFO=<ID=HM[>,]", header, re.M)),
        "VCF_CONTIG_CHR_PREFIX": bool(re.search(r"^##contig=<ID=chr", header, re.M)),
    }


def format_vcf_region(vcf_tags: dict[str, bool] | None, chrom: str,
                      start: int, end: int) -> str:
    norm = normalize_chromosome(chrom)
    prefix = "chr" if (vcf_tags and vcf_tags.get("VCF_CONTIG_CHR_PREFIX")) else ""
    return f"{prefix}{norm}:{start}-{end}"


def query_vcf(vcf: Path, region: str | None = None,
              minimum_lp: float | None = None,
              vcf_tags: dict[str, bool] | None = None) -> list[dict[str, Any]]:
    if vcf_tags is None:
        vcf_tags = inspect_vcf_tags(vcf)
    fmt = BASE_QUERY_FORMAT
    optional_cols: list[str] = []
    if vcf_tags.get("FORMAT_AF"):
        fmt += "\\t[%AF{0}]"
        optional_cols.append("AF")
    if vcf_tags.get("FORMAT_AFS"):
        fmt += "\\t[%AFS{0}]"
        optional_cols.append("AFS")
    if vcf_tags.get("INFO_HM"):
        fmt += "\\t%INFO/HM"
        optional_cols.append("HM")
    fmt += "\\n"

    if region is None and minimum_lp is not None:
        p_view = subprocess.Popen(
            ["bcftools", "view", "--threads", "4", "-u", "-i",
             f"FMT/LP[0:0]>={minimum_lp:.12g}", str(vcf)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        p_query = subprocess.Popen(
            ["bcftools", "query", "-f", fmt],
            stdin=p_view.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert p_view.stdout is not None
        p_view.stdout.close()
        stdout_text, q_err = p_query.communicate()
        v_err = p_view.stderr.read().decode("utf-8", errors="ignore") if p_view.stderr else ""
        if p_view.wait() != 0:
            raise subprocess.CalledProcessError(p_view.returncode, "bcftools view", stderr=v_err)
        if p_query.returncode != 0:
            raise subprocess.CalledProcessError(p_query.returncode, "bcftools query", stderr=q_err)
    else:
        command = ["bcftools", "query"]
        if region:
            command.extend(["-r", region])
        if minimum_lp is not None:
            command.extend(["-i", f"FMT/LP[0:0]>={minimum_lp:.12g}"])
        command.extend(["-f", fmt, str(vcf)])
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        stdout_text = result.stdout

    expected_len = 9 + len(optional_cols)
    records = []
    for line_number, line in enumerate(stdout_text.splitlines(), start=1):
        fields = line.split("\t")
        if len(fields) != expected_len:
            raise ValueError(f"Unexpected bcftools query row {line_number}: {line}")
        chromosome, position, rsid, ref, alt, lp, effect, se, sample_size = fields[:9]
        extra = dict(zip(optional_cols, fields[9:]))
        pos_int = int(position)
        chrom_norm = normalize_chromosome(chromosome)
        records.append({
            "CHROM": chrom_norm,
            "POS": pos_int,
            "ID": rsid,
            "REF": ref,
            "ALT": alt,
            "LP": numeric(lp),
            "ES": numeric(effect),
            "SE": numeric(se),
            "SS": numeric(sample_size),
            "AF": numeric(extra.get("AF")),
            "AFS": clean_str(extra.get("AFS")),
            "HM": clean_str(extra.get("HM")),
            "mhc_region": is_mhc_grch38(chrom_norm, pos_int),
            "mhc_core_grc": is_mhc_core_grc_grch38(chrom_norm, pos_int),
        })
    return records


# Backwards-compatible positional clustering helpers retained for unit tests
def cluster_significant(records: list[dict[str, Any]], merge_distance_bp: int
                        ) -> list[list[dict[str, Any]]]:
    significant = [record for record in records if record.get("LP") is not None]
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


def split_overwide_interval(
    chrom: str,
    start: int,
    end: int,
    leads: list[dict[str, Any]],
    sig_records: list[dict[str, Any]],
    max_width_bp: int,
) -> list[dict[str, Any]]:
    """Split a rare merged non-MHC interval wider than max_width_bp (<10 Mb Pan-UKB LD band)
    at the widest gap / lowest-density valley between greedy leads."""
    if end - start <= max_width_bp:
        best_lead = max(sig_records, key=lambda r: (r["LP"], -r["POS"]))
        return [{
            "CHROM": chrom,
            "start": start,
            "end": end,
            "leads": leads,
            "significant_records": sig_records,
            "lead": best_lead,
            "split_at_panukb_radius": False,
        }]

    sorted_leads = sorted(leads, key=lambda r: r["POS"])
    if len(sorted_leads) >= 2:
        # Choose cut point near the middle that maximizes inter-lead gap
        mid = (start + end) / 2.0
        best_cut = None
        best_score = -1.0
        for a, b in zip(sorted_leads, sorted_leads[1:]):
            if b["POS"] <= a["POS"]:
                continue
            cut = (a["POS"] + b["POS"]) // 2
            if cut <= start or cut >= end:
                continue
            gap = b["POS"] - a["POS"]
            balance = 1.0 - abs(cut - mid) / max(1.0, (end - start) / 2.0)
            score = gap * (0.5 + 0.5 * max(0.0, balance))
            if score > best_score:
                best_score = score
                best_cut = cut
        if best_cut is not None:
            left_sigs = [r for r in sig_records if r["POS"] <= best_cut]
            right_sigs = [r for r in sig_records if r["POS"] > best_cut]
            left_leads = [r for r in sorted_leads if r["POS"] <= best_cut] or left_sigs[:1]
            right_leads = [r for r in sorted_leads if r["POS"] > best_cut] or right_sigs[:1]
            if left_sigs and right_sigs:
                left_parts = split_overwide_interval(
                    chrom, start, best_cut, left_leads, left_sigs, max_width_bp
                )
                right_parts = split_overwide_interval(
                    chrom, best_cut + 1, end, right_leads, right_sigs, max_width_bp
                )
                for item in left_parts + right_parts:
                    item["split_at_panukb_radius"] = True
                return left_parts + right_parts

    # Fallback deterministic chunking if all significant SNPs sit in one tight spot
    cut = start + max_width_bp
    left_sigs = [r for r in sig_records if r["POS"] <= cut] or sig_records
    right_sigs = [r for r in sig_records if r["POS"] > cut]
    out = [{
        "CHROM": chrom,
        "start": start,
        "end": min(end, cut),
        "leads": sorted_leads,
        "significant_records": left_sigs,
        "lead": max(left_sigs, key=lambda r: (r["LP"], -r["POS"])),
        "split_at_panukb_radius": True,
    }]
    if right_sigs:
        out.extend(split_overwide_interval(
            chrom, cut + 1, end, right_sigs[:1], right_sigs, max_width_bp
        ))
    return out


def greedy_kanai_windows(
    records: list[dict[str, Any]],
    flank_bp: int = 1_500_000,
    max_region_width_bp: int = 9_500_000,
    mhc_bounds: tuple[int, int] = (MHC_START_GRCH38, MHC_END_GRCH38),
) -> list[dict[str, Any]]:
    """Kanai et al. (2021) greedy lead-centered window construction + overlap merge.

    1. Rank significant variants per chromosome by (-LP, POS).
    2. Select top unmasked variant, form [POS - flank_bp, POS + flank_bp] (clamped at MHC
       boundaries on chr6 so non-MHC flanks never bleed into the single MHC locus), and
       mask all significant variants inside that window from further lead selection.
    3. Merge overlapping windows and split any rare >max_region_width_bp chain so every
       locus stays strictly within Pan-UKB's 10 Mb stored LD radius.
    """
    if flank_bp < 1 or max_region_width_bp < 2 * flank_bp:
        raise ValueError("flank_bp must be >= 1 and max_region_width_bp >= 2 * flank_bp.")

    mhc_start, mhc_end = mhc_bounds
    by_chrom: dict[str, list[dict[str, Any]]] = {}
    for row in records:
        if row.get("LP") is None:
            continue
        chrom = normalize_chromosome(row["CHROM"])
        if chrom == MHC_CHROM and mhc_start <= row["POS"] <= mhc_end:
            continue
        by_chrom.setdefault(chrom, []).append(row)

    merged_regions: list[dict[str, Any]] = []
    for chrom in sorted(by_chrom, key=chromosome_sort_key):
        chrom_rows = sorted(by_chrom[chrom], key=lambda r: (-r["LP"], r["POS"]))
        unmasked = list(chrom_rows)
        raw_windows: list[dict[str, Any]] = []
        while unmasked:
            lead = unmasked[0]
            w_start = max(1, lead["POS"] - flank_bp)
            w_end = lead["POS"] + flank_bp
            if chrom == MHC_CHROM:
                if lead["POS"] < mhc_start:
                    w_end = min(w_end, mhc_start - 1)
                elif lead["POS"] > mhc_end:
                    w_start = max(w_start, mhc_end + 1)
            in_win = [r for r in unmasked if w_start <= r["POS"] <= w_end]
            unmasked = [r for r in unmasked if not (w_start <= r["POS"] <= w_end)]
            raw_windows.append({
                "start": w_start,
                "end": w_end,
                "leads": [lead],
                "sigs": in_win,
            })

        raw_windows.sort(key=lambda w: (w["start"], w["end"]))
        merged: list[dict[str, Any]] = []
        for win in raw_windows:
            if not merged or win["start"] > merged[-1]["end"]:
                merged.append({
                    "start": win["start"],
                    "end": win["end"],
                    "leads": list(win["leads"]),
                    "sigs": list(win["sigs"]),
                })
            else:
                merged[-1]["end"] = max(merged[-1]["end"], win["end"])
                merged[-1]["leads"].extend(win["leads"])
                merged[-1]["sigs"].extend(win["sigs"])

        for interval in merged:
            all_sigs = [
                r for r in chrom_rows
                if interval["start"] <= r["POS"] <= interval["end"]
            ]
            merged_regions.extend(split_overwide_interval(
                chrom=chrom,
                start=interval["start"],
                end=interval["end"],
                leads=interval["leads"],
                sig_records=all_sigs,
                max_width_bp=max_region_width_bp,
            ))

    merged_regions.sort(
        key=lambda r: (-r["lead"]["LP"], chromosome_sort_key(r["CHROM"]), r["start"])
    )
    return merged_regions


def build_analysis_regions(
    vcf: Path,
    vcf_tags: dict[str, bool],
    significant: list[dict[str, Any]],
    lifter: Any,
    significance_lp: float,
    cis_significance_lp: float | None,
    flank_bp: int,
    gene_window_bp: int,
    max_region_width_bp: int,
    exclude_mhc: bool,
    gene: dict[str, Any] | None,
    mhc_bounds: tuple[int, int] = (MHC_START_GRCH38, MHC_END_GRCH38),
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    """Construct the deduplicated list of regions to analyze (`cis_gene`, `mhc`, `lead_locus`)
    and structured QC summaries for `cis_gene` and `mhc`."""
    mhc_start, mhc_end = mhc_bounds
    cis_gate_lp = cis_significance_lp if cis_significance_lp is not None else significance_lp

    mhc_sigs = [
        r for r in significant
        if normalize_chromosome(r["CHROM"]) == MHC_CHROM and mhc_start <= r["POS"] <= mhc_end
    ]
    non_mhc_sigs = [
        r for r in significant
        if not (normalize_chromosome(r["CHROM"]) == MHC_CHROM and mhc_start <= r["POS"] <= mhc_end)
    ]

    cis_summary: dict[str, Any] = {"status": "not_applicable"}
    cis_region_spec: dict[str, Any] | None = None
    gene_overlaps_mhc = False

    if gene is not None:
        g_chrom = normalize_chromosome(gene["chrom"])
        cis_start = max(1, gene["start"] - gene_window_bp)
        cis_end = gene["end"] + gene_window_bp
        gene_overlaps_mhc = (
            g_chrom == MHC_CHROM
            and cis_start <= mhc_end
            and cis_end >= mhc_start
        )
        if gene_overlaps_mhc:
            cis_in_mhc = query_vcf(
                vcf,
                region=format_vcf_region(vcf_tags, MHC_CHROM, mhc_start, mhc_end),
                minimum_lp=cis_gate_lp,
                vcf_tags=vcf_tags,
            )
            max_cis_lp = max((r["LP"] for r in cis_in_mhc if r["LP"] is not None), default=None)
            if exclude_mhc:
                cis_summary = {
                    "status": "excluded_mhc_by_option",
                    "gene_symbol": gene["gene_symbol"],
                    "uniprot_id": gene["uniprot_id"],
                    "region_grch38": f"{g_chrom}:{cis_start}-{cis_end}",
                    "threshold_lp": cis_gate_lp,
                    "significant_variants_in_mhc": len(cis_in_mhc),
                    "max_lp": max_cis_lp,
                }
            elif cis_in_mhc:
                mhc_sigs = cis_in_mhc if len(cis_in_mhc) > len(mhc_sigs) else mhc_sigs
                cis_summary = {
                    "status": "merged_into_single_mhc_region",
                    "gene_symbol": gene["gene_symbol"],
                    "uniprot_id": gene["uniprot_id"],
                    "region_grch38": f"{MHC_CHROM}:{mhc_start}-{mhc_end}",
                    "threshold_lp": cis_gate_lp,
                    "significant_variants": len(cis_in_mhc),
                    "max_lp": max_cis_lp,
                }
            else:
                cis_summary = {
                    "status": "skipped_below_significance_threshold",
                    "gene_symbol": gene["gene_symbol"],
                    "uniprot_id": gene["uniprot_id"],
                    "region_grch38": f"{g_chrom}:{cis_start}-{cis_end}",
                    "threshold_lp": cis_gate_lp,
                    "significant_variants": 0,
                    "max_lp": None,
                }
        else:
            if g_chrom == MHC_CHROM:
                if gene["end"] < mhc_start:
                    cis_end = min(cis_end, mhc_start - 1)
                elif gene["start"] > mhc_end:
                    cis_start = max(cis_start, mhc_end + 1)
            cis_hits = (
                [r for r in non_mhc_sigs if r["CHROM"] == g_chrom and cis_start <= r["POS"] <= cis_end]
                if abs(cis_gate_lp - significance_lp) < 1e-9
                else query_vcf(
                    vcf,
                    region=format_vcf_region(vcf_tags, g_chrom, cis_start, cis_end),
                    minimum_lp=cis_gate_lp,
                    vcf_tags=vcf_tags,
                )
            )
            if cis_hits:
                lead = max(cis_hits, key=lambda r: (r["LP"], -r["POS"]))
                cis_region_spec = {
                    "region_id": f"cis_{gene['gene_symbol']}_chr{g_chrom}_{cis_start}_{cis_end}",
                    "region_type": "cis_gene",
                    "CHROM": g_chrom,
                    "start": cis_start,
                    "end": cis_end,
                    "lead": lead,
                    "significant_variant_count": len(cis_hits),
                    "mhc_region": False,
                    "mhc_core_grc": False,
                    "cis_overlap": True,
                    "split_at_panukb_radius": False,
                    "extended_by_adjacent_lead": False,
                }
                cis_summary = {
                    "status": "activated",
                    "gene_symbol": gene["gene_symbol"],
                    "uniprot_id": gene["uniprot_id"],
                    "region_id": cis_region_spec["region_id"],
                    "region_grch38": f"{g_chrom}:{cis_start}-{cis_end}",
                    "threshold_lp": cis_gate_lp,
                    "significant_variants": len(cis_hits),
                    "max_lp": lead["LP"],
                    "lead_id": lead["ID"],
                }
                # Mask genome-wide significant variants already inside cis_gene
                non_mhc_sigs = [
                    r for r in non_mhc_sigs
                    if not (r["CHROM"] == g_chrom and cis_start <= r["POS"] <= cis_end)
                ]
            else:
                # Query top LP in cis window lightly for QC audit
                cis_all = query_vcf(
                    vcf,
                    region=format_vcf_region(vcf_tags, g_chrom, cis_start, cis_end),
                    vcf_tags=vcf_tags,
                )
                max_cis_lp = max((r["LP"] for r in cis_all if r["LP"] is not None), default=None)
                cis_summary = {
                    "status": "skipped_below_significance_threshold",
                    "gene_symbol": gene["gene_symbol"],
                    "uniprot_id": gene["uniprot_id"],
                    "region_grch38": f"{g_chrom}:{cis_start}-{cis_end}",
                    "variants_in_window": len(cis_all),
                    "threshold_lp": cis_gate_lp,
                    "significant_variants": 0,
                    "max_lp": max_cis_lp,
                }

    # Build Kanai greedy 3 Mb merged windows for remaining non-MHC significant variants
    kanai_loci = greedy_kanai_windows(
        non_mhc_sigs,
        flank_bp=flank_bp,
        max_region_width_bp=max_region_width_bp,
        mhc_bounds=mhc_bounds,
    )

    regions: list[dict[str, Any]] = []
    # Merge any adjacent Kanai window that overlaps an active cis_gene (within max_region_width_bp)
    for kl in kanai_loci:
        if (
            cis_region_spec is not None
            and kl["CHROM"] == cis_region_spec["CHROM"]
            and kl["start"] <= cis_region_spec["end"]
            and kl["end"] >= cis_region_spec["start"]
        ):
            new_start = min(cis_region_spec["start"], kl["start"])
            new_end = max(cis_region_spec["end"], kl["end"])
            if new_end - new_start <= max_region_width_bp:
                cis_region_spec["start"] = new_start
                cis_region_spec["end"] = new_end
                cis_region_spec["region_id"] = (
                    f"cis_{gene['gene_symbol']}_chr{kl['CHROM']}_{new_start}_{new_end}"
                )
                cis_region_spec["significant_variant_count"] += len(kl["significant_records"])
                if kl["lead"]["LP"] > cis_region_spec["lead"]["LP"]:
                    cis_region_spec["lead"] = kl["lead"]
                cis_region_spec["extended_by_adjacent_lead"] = True
                cis_summary["region_id"] = cis_region_spec["region_id"]
                cis_summary["region_grch38"] = f"{kl['CHROM']}:{new_start}-{new_end}"
                cis_summary["extended_by_adjacent_lead"] = True
                continue
            # Otherwise trim non-cis window outside cis_gene so variants are not duplicated
            if kl["lead"]["POS"] < cis_region_spec["start"]:
                kl["end"] = cis_region_spec["start"] - 1
            else:
                kl["start"] = cis_region_spec["end"] + 1
            if kl["end"] <= kl["start"]:
                continue

        regions.append({
            "region_id": f"chr{kl['CHROM']}_{kl['start']}_{kl['end']}",
            "region_type": "lead_locus",
            "CHROM": kl["CHROM"],
            "start": kl["start"],
            "end": kl["end"],
            "lead": kl["lead"],
            "significant_variant_count": len(kl["significant_records"]),
            "mhc_region": False,
            "mhc_core_grc": False,
            "cis_overlap": False,
            "split_at_panukb_radius": kl["split_at_panukb_radius"],
            "extended_by_adjacent_lead": False,
        })

    if cis_region_spec is not None:
        regions.insert(0, cis_region_spec)

    # MHC handling
    mhc_summary: dict[str, Any]
    if mhc_sigs:
        mhc_lead = max(mhc_sigs, key=lambda r: (r["LP"], -r["POS"]))
        has_core = any(is_mhc_core_grc_grch38(r["CHROM"], r["POS"]) for r in mhc_sigs)
        if exclude_mhc:
            mhc_summary = {
                "status": "excluded_by_option",
                "region_grch38": f"{MHC_CHROM}:{mhc_start}-{mhc_end}",
                "core_grc_grch38": f"{MHC_CHROM}:{MHC_CORE_GRC_GRCH38[0]}-{MHC_CORE_GRC_GRCH38[1]}",
                "significant_variant_count": len(mhc_sigs),
                "overlaps_core_grc_mhc": has_core,
                "lead_id": mhc_lead["ID"],
                "lead_pos": mhc_lead["POS"],
                "lead_lp": mhc_lead["LP"],
            }
        else:
            mhc_region_id = f"chr{MHC_CHROM}_{mhc_start}_{mhc_end}_mhc"
            regions.append({
                "region_id": mhc_region_id,
                "region_type": "mhc",
                "CHROM": MHC_CHROM,
                "start": mhc_start,
                "end": mhc_end,
                "lead": mhc_lead,
                "significant_variant_count": len(mhc_sigs),
                "mhc_region": True,
                "mhc_core_grc": has_core,
                "cis_overlap": gene_overlaps_mhc,
                "split_at_panukb_radius": False,
                "extended_by_adjacent_lead": False,
            })
            mhc_summary = {
                "status": "included_single_region",
                "region_id": mhc_region_id,
                "region_grch38": f"{MHC_CHROM}:{mhc_start}-{mhc_end}",
                "core_grc_grch38": f"{MHC_CHROM}:{MHC_CORE_GRC_GRCH38[0]}-{MHC_CORE_GRC_GRCH38[1]}",
                "significant_variant_count": len(mhc_sigs),
                "overlaps_core_grc_mhc": has_core,
                "lead_id": mhc_lead["ID"],
                "lead_pos": mhc_lead["POS"],
                "lead_lp": mhc_lead["LP"],
            }
    else:
        mhc_summary = {
            "status": "no_significant_variants",
            "region_grch38": f"{MHC_CHROM}:{mhc_start}-{mhc_end}",
            "significant_variant_count": 0,
        }

    # Annotate lead liftover status on every candidate region
    for reg in regions:
        mapped, reason = map_variant(lifter, reg["lead"])
        reg["lead_GRCh37_CHROM"] = mapped["LD_CHROM"] if mapped else ""
        reg["lead_GRCh37_POS"] = mapped["LD_POS"] if mapped else ""
        reg["lead_liftover_status"] = "mapped" if mapped else "unmapped_lead"
        reg["lead_liftover_reason"] = reason or ""

    # Order: cis_gene first (if active), then by descending lead LP
    regions.sort(
        key=lambda r: (
            0 if r["region_type"] == "cis_gene" else 1,
            -r["lead"]["LP"],
            chromosome_sort_key(r["CHROM"]),
            r["start"],
        )
    )
    return regions, cis_summary, mhc_summary


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


def read_gene_annotation(path: Path, identifier: str,
                         strict: bool = True) -> dict[str, Any] | None:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    matches = [
        row for row in rows
        if row.get("uniprot_id") == identifier or row.get("gene_symbol") == identifier
    ]
    if len(matches) != 1:
        if strict:
            raise ValueError(
                f"Expected exactly one gene annotation for {identifier}; found {len(matches)}"
            )
        return None
    gene = dict(matches[0])
    if gene["genome_build"] != "GRCh38":
        raise ValueError("Gene annotation must use GRCh38 coordinates.")
    gene["chrom"] = normalize_chromosome(gene["chrom"])
    gene["start"] = int(gene["start"])
    gene["end"] = int(gene["end"])
    if gene["start"] < 1 or gene["end"] < gene["start"]:
        raise ValueError(f"Invalid gene coordinates for {identifier}")
    return gene


def resolve_target_gene(
    gene_path: Path | None,
    gwas_id: str,
    metadata: dict[str, Any],
    explicit_gene_id: str | None = None,
) -> tuple[dict[str, Any] | None, str]:
    if gene_path is None or not gene_path.is_file():
        if explicit_gene_id:
            raise FileNotFoundError(
                f"--gene-id {explicit_gene_id} requested, but --gene-annotation was not provided."
            )
        return None, "no_gene_annotation_supplied"

    if explicit_gene_id:
        gene = read_gene_annotation(gene_path, explicit_gene_id, strict=True)
        return gene, "resolved_explicit"

    candidates: list[str] = []
    for key in ("uniprot_id", "gene_symbol", "gene_id"):
        val = metadata.get(key)
        if isinstance(val, str) and val.strip():
            candidates.append(val.strip())
    token_match = re.match(r"^([A-Z0-9][A-Z0-9-]+?)(?:_|$)", gwas_id)
    if token_match:
        candidates.append(token_match.group(1))

    ambiguous_path = gene_path.with_name("gene_regions_GRCh38.ambiguous.tsv")
    ambiguous_ids: set[str] = set()
    if ambiguous_path.is_file():
        with ambiguous_path.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                if row.get("uniprot_id"):
                    ambiguous_ids.add(row["uniprot_id"])

    for cand in candidates:
        if cand in ambiguous_ids:
            return None, f"ambiguous_multi_gene_uniprot:{cand}"
        gene = read_gene_annotation(gene_path, cand, strict=False)
        if gene is not None:
            return gene, "resolved_from_metadata_or_id"

    return None, "non_gene_trait"


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


def r_runtime_environment(project_root: Path) -> dict[str, str]:
    environment = os.environ.copy()
    environment.setdefault("RENV_CONFIG_SYNCHRONIZED_CHECK", "FALSE")
    environment.setdefault("R_MAX_VSIZE", "64Gb")
    if environment.get("GWASFM_IN_CONTAINER") == "1":
        return environment
    tools_root = Path(
        environment.get(
            "GWASFM_TOOLS_ROOT", str(Path.home() / ".cache" / "gwasfinemapping")
        )
    ).expanduser()
    if not tools_root.is_absolute():
        tools_root = project_root / tools_root
    environment.setdefault("RENV_PROJECT", str(project_root))
    rprofile = project_root / ".Rprofile"
    if rprofile.is_file():
        environment.setdefault("R_PROFILE_USER", str(rprofile))
    if "RENV_PATHS_LIBRARY" not in environment:
        environment["RENV_PATHS_LIBRARY"] = str(tools_root / "renv-library")
    env_bin = tools_root / "env" / "bin"
    if env_bin.is_dir():
        environment["PATH"] = f"{env_bin}:{environment.get('PATH', '')}"
    return environment


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
        values = [record.get(field) for field in ("ES", "SE", "LP", "SS")]
        if any(value is None for value in values[:3]) or (values[1] is not None and values[1] <= 0):
            excluded.append({
                **record,
                "status": "excluded",
                "reason": "missing_or_invalid_ES_SE_or_LP",
            })
            continue
        lifted_by_id[unique_id] = lifted
        eligible.append((unique_id, record, lifted, values))

    output_path.parent.mkdir(parents=True, exist_ok=True)
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


def _load_fai_index(fasta: Path) -> dict[str, tuple[int, int, int, int]]:
    fai = Path(f"{fasta}.fai")
    if not fai.is_file():
        return {}
    index: dict[str, tuple[int, int, int, int]] = {}
    for line in fai.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 5:
            index[parts[0]] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))
    return index


def _fetch_fasta_bases(fh: Any, fai_entry: tuple[int, int, int, int],
                       pos_1based: int, length: int) -> str:
    seq_len, offset, line_bases, line_width = fai_entry
    if pos_1based < 1 or pos_1based + length - 1 > seq_len:
        return ""
    zero_pos = pos_1based - 1
    out: list[str] = []
    remaining = length
    cur = zero_pos
    while remaining > 0:
        line_idx, col_idx = divmod(cur, line_bases)
        take = min(remaining, line_bases - col_idx)
        byte_pos = offset + line_idx * line_width + col_idx
        fh.seek(byte_pos)
        out.append(fh.read(take).decode("ascii", errors="ignore").upper())
        cur += take
        remaining -= take
    return "".join(out)


def normalize_lifted_vcf(vcf: Path, fasta: Path,
                         output_path: Path | None = None
                         ) -> list[dict[str, Any]]:
    fai_map = _load_fai_index(fasta)
    pre_swapped_ids: set[str] = set()
    with tempfile.TemporaryDirectory(prefix="gwasfm-normalize-") as directory:
        norm_input = vcf
        if fai_map:
            adjusted_vcf = Path(directory) / "pre_norm.vcf"
            with fasta.open("rb") as fasta_fh, vcf.open("r", encoding="utf-8") as src, adjusted_vcf.open("w", encoding="utf-8") as dst:
                for line in src:
                    if line.startswith("#"):
                        dst.write(line)
                        continue
                    cols = line.rstrip("\n").split("\t")
                    chrom, pos_s, uid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
                    if len(ref) != len(alt) and chrom in fai_map:
                        pos_i = int(pos_s)
                        max_l = max(len(ref), len(alt))
                        ref_seq = _fetch_fasta_bases(fasta_fh, fai_map[chrom], pos_i, max_l)
                        if ref_seq and not ref_seq.startswith(ref) and ref_seq.startswith(alt):
                            cols[3], cols[4] = alt, ref
                            pre_swapped_ids.add(uid)
                            line = "\t".join(cols) + "\n"
                    dst.write(line)
            norm_input = adjusted_vcf

        normalized = Path(directory) / "normalized.vcf"
        run_checked([
            "bcftools", "norm", "-f", str(fasta), "-c", "w", "-Ov",
            "-o", str(normalized), str(norm_input),
        ])
        result = run_checked([
            "bcftools", "query", "-f",
            "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t"
            "[%ES{0}]\\t[%SE{0}]\\t[%LP{0}]\\t[%SS{0}]\\n",
            str(normalized),
        ])
    rows = []
    fai_checked = bool(fai_map)
    for line in result.stdout.splitlines():
        chrom, pos, unique_id, ref, alt, es, se, lp, ss = line.split("\t")
        was_pre_swapped = unique_id in pre_swapped_ids
        if was_pre_swapped:
            ref, alt = alt, ref
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
            "PRE_SWAPPED": was_pre_swapped,
            "FAI_CHECKED": fai_checked,
        })
    if output_path is not None:
        write_tsv(
            output_path,
            rows,
            ["source_variant_key", "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT",
             "ES", "SE", "LP", "SS"],
        )
    return rows


def match_panel_variants(index_path: Path, normalized: list[dict[str, Any]]) -> Any:
    """Vectorized VariantIndex lookup; prevents non-swapped tandem-repeat indels from cross-matching."""
    import numpy as np
    import pandas as pd
    from ldcov.io.variant_index import VariantIndex

    n = len(normalized)
    if n == 0:
        return pd.DataFrame({"idx": np.empty(0, dtype=np.int64), "flip": np.empty(0, dtype=bool)})

    vi = VariantIndex(str(index_path))
    candidates = pd.DataFrame({
        "contig": [row["LD_CHROM"] for row in normalized],
        "position": [row["LD_POS"] for row in normalized],
        "ref": [row["LD_REF"] for row in normalized],
        "alt": [row["LD_ALT"] for row in normalized],
    })
    swap_mask = np.array([
        len(row["LD_REF"]) == len(row["LD_ALT"])
        or bool(row.get("PRE_SWAPPED", False))
        or not bool(row.get("FAI_CHECKED", False))
        for row in normalized
    ], dtype=bool)

    idx_out = np.full(n, -1, dtype=np.int64)
    flip_out = np.zeros(n, dtype=bool)
    if np.any(swap_mask):
        res_swap = vi.match_variants(candidates.loc[swap_mask], allow_allele_swap=True)
        idx_out[swap_mask] = res_swap["idx"].to_numpy()
        flip_out[swap_mask] = res_swap["flip"].to_numpy()
    if np.any(~swap_mask):
        res_exact = vi.match_variants(candidates.loc[~swap_mask], allow_allele_swap=False)
        idx_out[~swap_mask] = res_exact["idx"].to_numpy()
        flip_out[~swap_mask] = res_exact["flip"].to_numpy()
    return pd.DataFrame({"idx": idx_out, "flip": flip_out})


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
    root = os.environ.get("GWASFM_REFERENCE_ROOT") or "~/.cache/gwasfinemapping/reference"
    return Path(root).expanduser() / "panukb"


def installed_population_index(population: str, panukb_dir: Path | None) -> Path | None:
    if panukb_dir is None:
        return None
    index = panukb_dir / f"UKBB.{population}.variant_index.parquet"
    sidecar = Path(f"{index}.json")
    if not index.is_file() or not sidecar.is_file():
        return None
    record = json.loads(sidecar.read_text(encoding="utf-8"))
    if (record.get("population"), record.get("region"), record.get("build")) != (population, "genome", "GRCh37"):
        raise ValueError(f"{sidecar} does not describe the whole-genome {population} index on GRCh37")
    return index


def cut_regional_index(whole: Path, chromosome: str, start: int, end: int, destination: Path) -> int:
    import pyarrow.parquet as pq

    table = pq.read_table(
        whole,
        filters=[
            ("contig", "=", chromosome),
            ("position", ">=", start),
            ("position", "<=", end),
        ],
    )
    temporary = destination.with_name(destination.name + ".tmp")
    pq.write_table(table, temporary)
    os.replace(temporary, destination)
    return table.num_rows


def ensure_regional_index(population: str, chromosome: str, start: int, end: int,
                          cache_dir: Path, scripts_dir: Path,
                          panukb_dir: Path | None = None) -> Path:
    region = f"{chromosome}:{start}-{end}"
    index = cache_dir / "panukb" / (
        f"UKBB.{population}.chr{chromosome}_{start}_{end}.variant_index.parquet"
    )
    sidecar = Path(f"{index}.json")
    expected = {"population": population, "region": region, "build": "GRCh37"}
    index.parent.mkdir(parents=True, exist_ok=True)
    if index.exists():
        if not sidecar.is_file() or json.loads(sidecar.read_text(encoding="utf-8")) != expected:
            raise ValueError(f"Existing variant index lacks matching provenance: {index}")
        return index
    whole = installed_population_index(population, panukb_dir)
    if whole is not None:
        rows = cut_regional_index(whole, chromosome, start, end, index)
        if rows == 0:
            index.unlink(missing_ok=True)
            raise ValueError(f"The installed {population} index has no variants in {region}")
    else:
        run_checked([
            "bash", str(scripts_dir / "build_panukb_variant_index.sh"),
            population, str(index), region,
        ])
    write_json(sidecar, expected)
    return index


def _normalize_ld_matrix_inplace(matrix: Any, ids: list[str]) -> dict[str, Any]:
    """Validate and normalize a symmetric float64 LD matrix in-place using O(block^2) extra RAM."""
    import numpy as np

    n = len(ids)
    if matrix.shape != (n, n):
        raise ValueError("Extracted LD dimensions do not match its variant table.")
    if len(set(ids)) != n:
        raise ValueError("Extracted LD contains duplicate variant identifiers.")
    if n == 0:
        return {
            "variants": 0,
            "shape": [0, 0],
            "raw_diagonal_min": 0.0,
            "raw_diagonal_max": 0.0,
            "diagonal_normalization": "D^-1/2 R D^-1/2",
            "max_abs_asymmetry": 0.0,
            "missing_cells": 0,
        }

    step = 512
    for i in range(0, n, step):
        if not bool(np.all(np.isfinite(matrix[i:min(n, i + step), :]))):
            raise ValueError("Extracted LD has missing values or inconsistent indices.")

    diagonal = np.diag(matrix).copy()
    if np.any(diagonal <= 0) or np.max(np.abs(diagonal - 1.0)) > 0.01:
        raise ValueError("Pan-UKB LD diagonal is outside the 0.01 normalization tolerance.")

    max_asym = 0.0
    for i in range(0, n, step):
        ie = min(n, i + step)
        for j in range(i, n, step):
            je = min(n, j + step)
            diff = float(np.max(np.abs(matrix[i:ie, j:je] - matrix[j:je, i:ie].T)))
            if diff > max_asym:
                max_asym = diff
    if max_asym > 1e-8:
        raise ValueError("Extracted LD matrix is not symmetric.")

    scale = np.sqrt(diagonal)
    matrix /= scale[:, None]
    matrix /= scale[None, :]

    max_abs = 0.0
    for i in range(0, n, step):
        ie = min(n, i + step)
        for j in range(i, n, step):
            je = min(n, j + step)
            if i == j:
                blk = matrix[i:ie, i:ie]
                avg = 0.5 * (blk + blk.T)
                matrix[i:ie, i:ie] = avg
            else:
                avg = 0.5 * (matrix[i:ie, j:je] + matrix[j:je, i:ie].T)
                matrix[i:ie, j:je] = avg
                matrix[j:je, i:ie] = avg.T
            blk_max = float(np.max(np.abs(matrix[i:ie, j:je])))
            if blk_max > max_abs:
                max_abs = blk_max
        if not bool(np.all(np.isfinite(matrix[i:ie, :]))):
            raise ValueError("Normalized LD contains invalid correlations.")

    np.fill_diagonal(matrix, 1.0)
    if max_abs > 1.00001:
        raise ValueError("Normalized LD contains invalid correlations.")

    return {
        "variants": n,
        "shape": [n, n],
        "raw_diagonal_min": float(diagonal.min()),
        "raw_diagonal_max": float(diagonal.max()),
        "diagonal_normalization": "D^-1/2 R D^-1/2",
        "max_abs_asymmetry": max_asym,
        "missing_cells": 0,
    }


def _write_symmetric_fortran_npy(matrix: Any, ids: list[str], output_npy_path: Path) -> None:
    """Write a validated symmetric C- or F-contiguous float64 matrix as Fortran-order .npy with zero copy."""
    import numpy as np

    output_npy_path.parent.mkdir(parents=True, exist_ok=True)
    if matrix.flags["F_CONTIGUOUS"]:
        view_f = matrix
    elif matrix.flags["C_CONTIGUOUS"] and matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1]:
        # Symmetric matrix: transpose view shares buffer with zero copy and is F_CONTIGUOUS!
        view_f = matrix.T
    else:
        view_f = np.asfortranarray(matrix, dtype="<f8")
    np.save(output_npy_path, view_f, allow_pickle=False)
    output_npy_path.with_suffix(".ids.txt").write_text("\n".join(ids) + "\n", encoding="utf-8")


def _validate_and_normalize_extracted_ld(
    npz_path: Path, variants_path: Path
) -> tuple[Any, list[str], Any, dict[str, Any]]:
    import numpy as np
    import pandas as pd

    ld_variants = pd.read_csv(variants_path, sep="\t", dtype={"chrom": str})
    with np.load(npz_path) as data:
        matrix = np.array(data["matrix"], dtype="<f8", copy=True)
        indexes = np.asarray(data["idx"])
    if len(indexes) != len(ld_variants):
        raise ValueError("Extracted LD has missing values or inconsistent indices.")
    ids = [
        f"{normalize_chromosome(row.chrom)}:{int(row.pos)}:{row.ref}:{row.alt}"
        for row in ld_variants.itertuples(index=False)
    ]
    summary = _normalize_ld_matrix_inplace(matrix, ids)
    view_f = matrix.T if (matrix.flags["C_CONTIGUOUS"] and matrix.shape[0] == matrix.shape[1]) else np.asfortranarray(matrix, dtype="<f8")
    return view_f, ids, indexes, summary


def matrix_to_npy(npz_path: Path, variants_path: Path, output_npy_path: Path
                  ) -> dict[str, Any]:
    """Validate and write binary float64 .npy LD matrix + companion .ids.txt for fast R load."""
    normalized, ids, _, summary = _validate_and_normalize_extracted_ld(
        npz_path, variants_path
    )
    _write_symmetric_fortran_npy(normalized, ids, output_npy_path)
    return summary


# Retained for backwards compatibility with existing unit tests
def matrix_to_tsv(npz_path: Path, variants_path: Path, output_path: Path
                  ) -> dict[str, Any]:
    import numpy as np

    normalized, ids, indexes, summary = _validate_and_normalize_extracted_ld(
        npz_path, variants_path
    )
    normalized_npz = output_path.with_name("normalized_matrix.npz")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(normalized_npz, matrix=normalized, idx=indexes)
    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["LD_VARIANT_ID", *ids])
        for variant_id, row in zip(ids, normalized):
            writer.writerow([variant_id, *(f"{value:.12g}" for value in row)])
    return summary


def extract_ld_from_blockmatrix(
    bm_path: str,
    matched_rows: list[dict[str, Any]],
    ld_prefix: Path,
    block_cache: int = 4,
    fetch_workers: int = 8,
    storage_options: dict[str, Any] | None = None,
    direct_npy_path: Path | None = None,
    keep_npz: bool = True,
) -> dict[str, Any] | None:
    """Extract submatrix from Pan-UKB Hail BlockMatrix using pre-resolved PANEL_IDX + PANEL_ALLELE_SWAPPED."""
    import gc
    import numpy as np
    import pandas as pd
    from ldcov.io.blockmatrix.reader import HailBlockMatrixReader
    from ldcov.ld_bm.extract import _apply_flips, _assemble

    global_idxs = [int(row["PANEL_IDX"]) for row in matched_rows]
    if len(set(global_idxs)) != len(global_idxs):
        raise ValueError("Duplicate PANEL_IDX in matched_rows before BlockMatrix extraction.")

    reader = HailBlockMatrixReader(
        bm_path, storage_options=storage_options, block_cache=block_cache
    )
    matrix = _assemble(reader, global_idxs, max_workers=fetch_workers)
    del reader
    flip_positions = [i for i, row in enumerate(matched_rows) if row["PANEL_ALLELE_SWAPPED"]]
    matrix = _apply_flips(matrix, flip_positions)

    if keep_npz or direct_npy_path is None:
        targets = pd.DataFrame({
            "idx": np.array(global_idxs, dtype=np.int64),
            "chrom": [str(row["LD_CHROM"]) for row in matched_rows],
            "pos": [int(row["LD_POS"]) for row in matched_rows],
            "ref": [str(row["LD_REF"]) for row in matched_rows],
            "alt": [str(row["LD_ALT"]) for row in matched_rows],
            "rsid": [f"v{i:08d}" for i in range(len(matched_rows))],
            "flipped": [bool(row["PANEL_ALLELE_SWAPPED"]) for row in matched_rows],
        })
        ld_prefix.parent.mkdir(parents=True, exist_ok=True)
        targets.to_csv(f"{ld_prefix}.variants.tsv", sep="\t", index=False)
        np.savez(f"{ld_prefix}.npz", matrix=matrix, idx=targets["idx"].to_numpy())

    if direct_npy_path is not None:
        ids = [
            f"{normalize_chromosome(str(row['LD_CHROM']))}:{int(row['LD_POS'])}:{row['LD_REF']}:{row['LD_ALT']}"
            for row in matched_rows
        ]
        summary = _normalize_ld_matrix_inplace(matrix, ids)
        _write_symmetric_fortran_npy(matrix, ids, direct_npy_path)
        del matrix
        gc.collect()
        return summary

    del matrix
    gc.collect()
    return None


def extract_ld_from_local_panel(
    local_panel_npz: Path,
    matched_rows: list[dict[str, Any]],
    ld_prefix: Path,
) -> None:
    """Offline fixture LD extractor mirroring ldcov output format for test.sh / benchmark.sh."""
    import numpy as np

    with np.load(local_panel_npz) as panel:
        panel_matrix = np.asarray(panel["matrix"], dtype=np.float64)
        panel_idx = [int(x) for x in panel["idx"]]
    pos_by_idx = {idx: pos for pos, idx in enumerate(panel_idx)}
    missing = [row["PANEL_IDX"] for row in matched_rows if row["PANEL_IDX"] not in pos_by_idx]
    if missing:
        raise ValueError(f"Offline LD panel {local_panel_npz} is missing panel idx: {missing[:5]}")

    positions = [pos_by_idx[row["PANEL_IDX"]] for row in matched_rows]
    sub = panel_matrix[np.ix_(positions, positions)].copy()
    signs = np.array([-1.0 if row["PANEL_ALLELE_SWAPPED"] else 1.0 for row in matched_rows])
    sub = sub * signs[:, None] * signs[None, :]

    ld_prefix.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        f"{ld_prefix}.npz",
        matrix=sub,
        idx=np.array([row["PANEL_IDX"] for row in matched_rows], dtype=np.int64),
    )
    with Path(f"{ld_prefix}.variants.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["idx", "chrom", "pos", "ref", "alt", "rsid", "flipped"])
        for i, row in enumerate(matched_rows):
            writer.writerow([
                row["PANEL_IDX"],
                row["LD_CHROM"],
                row["LD_POS"],
                row["LD_REF"],
                row["LD_ALT"],
                f"v{i:08d}",
                row["PANEL_ALLELE_SWAPPED"],
            ])


def run_methods(
    variants_path: Path,
    matrix_path: Path | None,
    output_dir: Path,
    reference_n: int | None,
    region_id: str,
    region_type: str,
    methods: str = METHODS,
    write_tsv_sidecars: bool = False,
    effective_n: float | None = None,
    max_causal: int = 10,
    max_iter: int = 200,
) -> dict[str, Any]:
    project_root = Path(__file__).resolve().parents[1]
    command = [
        "Rscript", str(Path(__file__).with_name("run_region_methods.R")),
        "--variants", str(variants_path.resolve()),
        "--outdir", str(output_dir.resolve()),
        "--methods", methods,
        "--region-id", region_id,
        "--region-type", region_type,
        "--write-tsv", "true" if write_tsv_sidecars else "false",
        "--max-causal", str(max_causal),
        "--max-iter", str(max_iter),
    ]
    if reference_n is not None:
        command.extend(["--reference-sample-size", str(reference_n)])
    if effective_n is not None:
        command.extend(["--effective-n", str(effective_n)])
    if matrix_path is not None:
        command.extend(["--matrix", str(matrix_path.resolve())])
    run_checked(command, env=r_runtime_environment(project_root), cwd=project_root)
    qc_file = output_dir / "method_qc.json"
    return json.loads(qc_file.read_text(encoding="utf-8")) if qc_file.is_file() else {}


def aggregate_method_results(
    output_root: Path,
    scripts_dir: Path,
    input_root: Path | None = None,
    allow_empty: bool = True,
) -> None:
    project_root = scripts_dir.parent
    command = [
        "Rscript", str(scripts_dir / "aggregate_method_results.R"),
        "--output-root", str(output_root.resolve()),
        "--allow-empty", "true" if allow_empty else "false",
        "--compress", "xz",
    ]
    if input_root is not None:
        command.extend(["--input-root", str(input_root.resolve())])
    run_checked(command, env=r_runtime_environment(project_root), cwd=project_root)


def render_gwas_report(output_dir: Path, scripts_dir: Path) -> None:
    report_script = scripts_dir / "render_report.R"
    if not report_script.is_file():
        return
    project_root = scripts_dir.parent
    command = [
        "Rscript", str(report_script),
        "--gwas-dir", str(output_dir.resolve()),
    ]
    run_checked(command, env=r_runtime_environment(project_root), cwd=project_root)


def analyze_locus(
    locus: dict[str, Any],
    vcf: Path,
    vcf_tags: dict[str, bool],
    lifter: Any,
    fasta: Path,
    population: str,
    reference_sample_size: int,
    cache_dir: Path,
    work_locus_dir: Path,
    rds_stage_dir: Path,
    scripts_dir: Path,
    gene: dict[str, Any] | None,
    gene_window_bp: int,
    panukb_dir: Path | None = None,
    local_ld_panel: Path | None = None,
    keep_intermediates: bool = False,
    effective_n: float | None = None,
    max_causal: int = 10,
    max_iter: int = 200,
) -> dict[str, Any]:
    t0 = time.monotonic()
    chrom = locus["CHROM"]
    start = locus["start"]
    end = locus["end"]
    region_id = locus["region_id"]
    region_type = locus["region_type"]
    region_str = format_vcf_region(vcf_tags, chrom, start, end)

    records = query_vcf(vcf, region=region_str, vcf_tags=vcf_tags)
    if not records:
        return {
            "region_id": region_id,
            "region_type": region_type,
            "source_region_grch38": f"{chrom}:{start}-{end}",
            "status": "skipped_empty_vcf_region",
            "mhc_region": locus["mhc_region"],
            "mhc_core_grc": locus["mhc_core_grc"],
            "variants_in_source_region": 0,
        }

    inputs_dir = work_locus_dir / "inputs"
    mapping_dir = work_locus_dir / "mapping"
    ld_dir = work_locus_dir / "ld"
    methods_dir = work_locus_dir / "methods"
    for d in (inputs_dir, mapping_dir, ld_dir, methods_dir):
        d.mkdir(parents=True, exist_ok=True)

    lifted_vcf = inputs_dir / "lifted.vcf"
    lifted_by_id, liftover_excluded = write_lifted_vcf(records, lifter, lifted_vcf)
    if not lifted_by_id:
        reason_counts = dict(Counter(r["reason"] for r in liftover_excluded))
        return {
            "region_id": region_id,
            "region_type": region_type,
            "source_region_grch38": f"{chrom}:{start}-{end}",
            "status": "skipped_no_variants_mapped_grch37",
            "mhc_region": locus["mhc_region"],
            "mhc_core_grc": locus["mhc_core_grc"],
            "variants_in_source_region": len(records),
            "variants_lifted": 0,
            "variants_ld_matched": 0,
            "exclusion_reasons": reason_counts,
        }

    normalized = normalize_lifted_vcf(
        lifted_vcf,
        fasta,
        mapping_dir / "normalized_grch37_variants.tsv" if keep_intermediates else None,
    )
    normalized = [
        {**lifted_by_id[row["source_variant_key"]], **row}
        for row in normalized
        if row["source_variant_key"] in lifted_by_id
    ]
    # Keep only variants mapping onto the expected target chromosome
    normalized_same_chrom = [r for r in normalized if r["LD_CHROM"] == chrom]
    cross_chrom_excluded = [
        {**r, "status": "excluded", "reason": "lifted_to_different_chromosome"}
        for r in normalized
        if r["LD_CHROM"] != chrom
    ]
    normalized = normalized_same_chrom
    if not normalized:
        reason_counts = dict(Counter(
            r["reason"] for r in (liftover_excluded + cross_chrom_excluded)
        ))
        return {
            "region_id": region_id,
            "region_type": region_type,
            "source_region_grch38": f"{chrom}:{start}-{end}",
            "status": "skipped_no_same_chrom_grch37_variants",
            "mhc_region": locus["mhc_region"],
            "mhc_core_grc": locus["mhc_core_grc"],
            "variants_in_source_region": len(records),
            "variants_lifted": 0,
            "variants_ld_matched": 0,
            "exclusion_reasons": reason_counts,
        }

    target_chrom = chrom
    target_start = max(1, min(row["LD_POS"] for row in normalized) - 100)
    target_end = max(row["LD_POS"] for row in normalized) + 100
    grch37_span = target_end - target_start
    if grch37_span > 9_990_000:
        raise ValueError(
            f"Region {region_id} spans {grch37_span:,} bp on GRCh37 ({target_chrom}:{target_start}-{target_end}), "
            "exceeding Pan-UKB's 10 Mb stored LD radius."
        )

    t_idx0 = time.monotonic()
    index = ensure_regional_index(
        population, target_chrom, target_start, target_end, cache_dir, scripts_dir, panukb_dir
    )
    index_seconds = round(time.monotonic() - t_idx0, 3)

    matches = match_panel_variants(index, normalized)
    matched, panel_excluded = resolve_panel_matches(normalized, matches)
    all_excluded = liftover_excluded + cross_chrom_excluded + panel_excluded
    reason_counts = dict(Counter(r["reason"] for r in all_excluded))

    if not matched:
        return {
            "region_id": region_id,
            "region_type": region_type,
            "source_region_grch38": f"{chrom}:{start}-{end}",
            "target_region_grch37": f"{target_chrom}:{target_start}-{target_end}",
            "status": "skipped_no_panel_matched_variants",
            "mhc_region": locus["mhc_region"],
            "mhc_core_grc": locus["mhc_core_grc"],
            "variants_in_source_region": len(records),
            "variants_lifted": len(normalized),
            "variants_ld_matched": 0,
            "exclusion_reasons": reason_counts,
        }

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
            cis_trans_label = (
                "cis"
                if gene is not None
                and source["CHROM"] == gene["chrom"]
                and (gene["start"] - gene_window_bp) <= source["POS"] <= (gene["end"] + gene_window_bp)
                else ("trans" if gene is not None else "na")
            )
            gwas_variant_id = (
                f"{source['CHROM']}:{source['POS']}:{source['REF']}:{source['ALT']}"
            )
            out = {
                **source,
                "GWAS_VARIANT_ID": gwas_variant_id,
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
                "cis_trans": cis_trans_label,
                "mhc_region": is_mhc_grch38(source["CHROM"], source["POS"]),
                "mhc_core_grc": (
                    is_mhc_core_grc_grch38(source["CHROM"], source["POS"])
                    or (
                        row["LD_CHROM"] == MHC_CHROM
                        and MHC_CORE_GRC_GRCH37[0] <= row["LD_POS"] <= MHC_CORE_GRC_GRCH37[1]
                    )
                ),
            }
            matched_rows.append(out)
            stream.write(
                f"v{i:08d} {row['LD_CHROM']} {row['LD_POS']} "
                f"{row['LD_REF']} {row['LD_ALT']}\n"
            )

    if keep_intermediates:
        mapping_rows = [
            {**row, "status": "matched", "reason": ""}
            for row in matched
        ] + all_excluded
        mapping_columns = [
            "CHROM", "POS", "ID", "REF", "ALT", "LP", "ES", "SE", "SS", "AF", "AFS", "HM",
            "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "LIFTOVER_STRAND",
            "LIFTOVER_CHAIN", "PANEL_IDX", "PANEL_ALLELE_SWAPPED", "status", "reason",
            "duplicate_of",
        ]
        write_tsv(mapping_dir / "variant_mapping.tsv", mapping_rows, mapping_columns)

    matched_variants_path = inputs_dir / "matched_variants.tsv"
    write_tsv(
        matched_variants_path,
        matched_rows,
        [
            "CHROM", "POS", "ID", "REF", "ALT",
            "GWAS_VARIANT_ID", "GWAS_CHROM", "GWAS_POS", "GWAS_REF", "GWAS_ALT",
            "LP", "ES", "SE", "SS", "AF", "AFS", "HM", "EZ",
            "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "LD_VARIANT_ID",
            "PANEL_IDX", "PANEL_ALLELE_SWAPPED", "cis_trans",
            "mhc_region", "mhc_core_grc",
        ],
    )

    ld_prefix = ld_dir / f"panukb_{population.lower()}"
    matrix_npy_path = ld_dir / "matrix.npy"
    if local_ld_panel is not None:
        extract_ld_from_local_panel(local_ld_panel, matched_rows, ld_prefix)
        ld_validation = matrix_to_npy(
            Path(f"{ld_prefix}.npz"),
            Path(f"{ld_prefix}.variants.tsv"),
            matrix_npy_path,
        )
        if not keep_intermediates:
            Path(f"{ld_prefix}.npz").unlink(missing_ok=True)
            Path(f"{ld_prefix}.variants.tsv").unlink(missing_ok=True)
    else:
        ld_val_opt = extract_ld_from_blockmatrix(
            bm_path=f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.bm",
            matched_rows=matched_rows,
            ld_prefix=ld_prefix,
            block_cache=4,
            fetch_workers=8,
            direct_npy_path=matrix_npy_path,
            keep_npz=keep_intermediates,
        )
        assert ld_val_opt is not None
        ld_validation = ld_val_opt
    if keep_intermediates:
        write_json(ld_dir / "validation.json", ld_validation)

    method_qc = run_methods(
        variants_path=matched_variants_path,
        matrix_path=matrix_npy_path,
        output_dir=methods_dir,
        reference_n=reference_sample_size,
        region_id=region_id,
        region_type=region_type,
        write_tsv_sidecars=keep_intermediates,
        effective_n=effective_n,
        max_causal=max_causal,
        max_iter=max_iter,
    )

    # Stage compact per-method region_result.rds files for final .rda aggregation
    for m_id in METHODS.split(","):
        src_rds = methods_dir / m_id / "region_result.rds"
        if src_rds.is_file():
            dst_rds = rds_stage_dir / region_id / m_id / "region_result.rds"
            dst_rds.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src_rds), str(dst_rds))

    # Delete regional cut index in cache if huge? Keep small parquet slice in cache_dir,
    # and clean locus scratch directory immediately unless --keep-intermediates is requested.
    if not keep_intermediates:
        shutil.rmtree(work_locus_dir, ignore_errors=True)

    n_swapped = sum(1 for r in matched_rows if r["PANEL_ALLELE_SWAPPED"])
    lead_matched = max(matched_rows, key=lambda r: (r["LP"], -r["GWAS_POS"]))
    ss_vals = sorted(float(r["SS"]) for r in matched_rows if r.get("SS") is not None)
    ss_min = ss_vals[0] if ss_vals else None
    ss_max = ss_vals[-1] if ss_vals else None
    ss_median = ss_vals[len(ss_vals) // 2] if ss_vals else None
    ss_ratio = round(ss_min / ss_max, 4) if (ss_min is not None and ss_max and ss_max > 0) else None
    ss_het_warn = bool(ss_ratio is not None and ss_ratio < 0.80)
    return {
        "region_id": region_id,
        "region_type": region_type,
        "status": "analyzed",
        "source_region_grch38": f"{chrom}:{start}-{end}",
        "target_region_grch37": f"{target_chrom}:{target_start}-{target_end}",
        "width_bp_grch38": end - start,
        "width_bp_grch37": grch37_span,
        "mhc_region": locus["mhc_region"],
        "mhc_core_grc": any(r["mhc_core_grc"] for r in matched_rows),
        "cis_overlap": locus["cis_overlap"],
        "split_at_panukb_radius": locus.get("split_at_panukb_radius", False),
        "extended_by_adjacent_lead": locus.get("extended_by_adjacent_lead", False),
        "significant_variants_in_region": locus["significant_variant_count"],
        "variants_in_source_region": len(records),
        "variants_lifted": len(normalized),
        "variants_ld_matched": len(matched_rows),
        "variants_allele_swapped": n_swapped,
        "variants_excluded": len(records) - len(matched_rows),
        "exclusion_reasons": reason_counts,
        "ss_min": ss_min,
        "ss_median": ss_median,
        "ss_max": ss_max,
        "ss_min_max_ratio": ss_ratio,
        "sample_size_heterogeneity_warning": ss_het_warn,
        "lead_variant_id": lead_matched["ID"],
        "lead_gwas_variant_id": lead_matched["GWAS_VARIANT_ID"],
        "lead_pos_grch38": lead_matched["GWAS_POS"],
        "lead_lp": lead_matched["LP"],
        "index_cut_seconds": index_seconds,
        "wall_seconds": round(time.monotonic() - t0, 3),
        "ld_validation": ld_validation,
        "methods": method_qc,
    }


def load_or_init_sidecars(
    vcf: Path,
    metadata_arg: str | None,
    qc_arg: str | None,
    study_id_arg: str | None,
    gwas_id_arg: str | None,
    trait_arg: str | None,
    ancestry_arg: str | None,
) -> tuple[dict[str, Any], dict[str, Any], str, str]:
    metadata: dict[str, Any] = {}
    qc_upstream: dict[str, Any] = {}
    if metadata_arg:
        mp = Path(metadata_arg).resolve()
        if not mp.is_file():
            raise FileNotFoundError(f"Metadata JSON not found: {mp}")
        metadata = json.loads(mp.read_text(encoding="utf-8"))
    if qc_arg:
        qp = Path(qc_arg).resolve()
        if not qp.is_file():
            raise FileNotFoundError(f"QC JSON not found: {qp}")
        qc_upstream = json.loads(qp.read_text(encoding="utf-8"))

    study_id = study_id_arg or metadata.get("study_id") or "study"
    default_gwas_id = re.sub(r"\.vcf(\.gz|\.bgz)?$", "", vcf.name)
    gwas_id = gwas_id_arg or metadata.get("output_id") or metadata.get("gwas_id") or default_gwas_id

    coord_meta = metadata.get("coordinate_build", "GRCh38")
    coord_qc = qc_upstream.get("coordinate_build", "GRCh38")
    if coord_meta != "GRCh38" or coord_qc != "GRCh38":
        raise ValueError("GWASfinemapping requires a GRCh38 GWAS-VCF.")

    metadata.setdefault("study_id", study_id)
    metadata.setdefault("output_id", gwas_id)
    metadata.setdefault("coordinate_build", "GRCh38")
    if trait_arg:
        metadata["phenotype"] = trait_arg
    metadata.setdefault("phenotype", gwas_id)
    if ancestry_arg:
        metadata["ancestry"] = ancestry_arg
    return metadata, qc_upstream, study_id, gwas_id


def run(args: argparse.Namespace) -> Path:
    t_run0 = time.monotonic()
    vcf = Path(args.vcf).resolve()
    vcf_tbi = Path(f"{vcf}.tbi")
    vcf_csi = Path(f"{vcf}.csi")
    chain = Path(args.chain).resolve()
    fasta = Path(args.grch37_fasta).resolve()
    fasta_index = Path(f"{fasta}.fai")
    gene_path = Path(args.gene_annotation).resolve() if args.gene_annotation else None
    local_ld_panel = Path(args.local_ld_panel).resolve() if args.local_ld_panel else None
    scripts_dir = Path(__file__).resolve().parent

    if not all(path.is_file() for path in (vcf, chain, fasta, fasta_index)):
        raise FileNotFoundError("VCF, chain, GRCh37 FASTA, and FASTA .fai must exist.")
    if not (vcf_tbi.is_file() or vcf_csi.is_file()):
        raise FileNotFoundError(f"Tabix/CSI index (.tbi or .csi) missing for {vcf}")

    metadata, qc_upstream, study_id, gwas_id = load_or_init_sidecars(
        vcf=vcf,
        metadata_arg=args.metadata,
        qc_arg=args.qc,
        study_id_arg=args.study_id,
        gwas_id_arg=args.gwas_id,
        trait_arg=args.trait,
        ancestry_arg=args.ancestry,
    )
    gene, gene_resolution = resolve_target_gene(
        gene_path=gene_path,
        gwas_id=gwas_id,
        metadata=metadata,
        explicit_gene_id=args.gene_id,
    )
    reference_n = resolve_reference_sample_size(args.population, args.reference_sample_size)
    if reference_n is None:
        raise ValueError(f"Supply --reference-sample-size for the {args.population} LD panel.")

    gwas_ancestry = metadata.get("ancestry")
    if (
        gwas_ancestry in REFERENCE_SAMPLES
        and gwas_ancestry != args.population
        and not getattr(args, "allow_ancestry_mismatch", False)
    ):
        raise ValueError(
            f"GWAS ancestry {gwas_ancestry} does not match LD population {args.population}."
        )
    if gwas_ancestry == "ALL":
        ancestry_ld_status = "multi_ancestry_gwas_single_panel_ld"
    elif gwas_ancestry in REFERENCE_SAMPLES and gwas_ancestry != args.population:
        ancestry_ld_status = "explicit_override_ancestry_mismatch"
    elif gwas_ancestry == args.population:
        ancestry_ld_status = "matched"
    else:
        ancestry_ld_status = "unspecified_gwas_ancestry"

    output_dir = Path(args.output_root).resolve() / study_id / gwas_id
    if output_dir.exists() and any(output_dir.iterdir()):
        if getattr(args, "force", False):
            shutil.rmtree(output_dir)
        else:
            raise FileExistsError(f"GWAS output directory is not empty: {output_dir} (pass --force to overwrite)")
    output_dir.mkdir(parents=True, exist_ok=True)

    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
    vcf_tags = inspect_vcf_tags(vcf)
    significant = query_vcf(vcf, minimum_lp=args.significance_lp, vcf_tags=vcf_tags)
    significant.sort(key=lambda r: (-r["LP"], chromosome_sort_key(r["CHROM"]), r["POS"]))

    candidate_regions, cis_summary, mhc_summary = build_analysis_regions(
        vcf=vcf,
        vcf_tags=vcf_tags,
        significant=significant,
        lifter=lifter,
        significance_lp=args.significance_lp,
        cis_significance_lp=args.cis_significance_lp,
        flank_bp=args.region_flank_bp,
        gene_window_bp=args.gene_window_bp,
        max_region_width_bp=args.max_region_width_bp,
        exclude_mhc=args.exclude_mhc,
        gene=gene,
    )

    selected_regions = (
        candidate_regions[: args.max_loci] if args.max_loci > 0 else candidate_regions
    )
    panukb_dir = (
        Path(args.panukb_dir).expanduser().resolve()
        if args.panukb_dir
        else default_panukb_dir()
    )
    cache_dir = Path(args.cache_dir).expanduser().resolve()

    analyzed: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix=f"gwasfm-{gwas_id[:24]}-") as scratch_tmp:
        scratch_root = Path(scratch_tmp)
        rds_stage_dir = scratch_root / "rds_stage"
        rds_stage_dir.mkdir(parents=True, exist_ok=True)

        for locus in selected_regions:
            work_locus_dir = (
                output_dir / "loci" / locus["region_id"]
                if args.keep_intermediates
                else scratch_root / "loci" / locus["region_id"]
            )
            locus_report = analyze_locus(
                locus=locus,
                vcf=vcf,
                vcf_tags=vcf_tags,
                lifter=lifter,
                fasta=fasta,
                population=args.population,
                reference_sample_size=reference_n,
                cache_dir=cache_dir,
                work_locus_dir=work_locus_dir,
                rds_stage_dir=rds_stage_dir,
                scripts_dir=scripts_dir,
                gene=gene,
                gene_window_bp=args.gene_window_bp,
                panukb_dir=panukb_dir,
                local_ld_panel=local_ld_panel,
                keep_intermediates=args.keep_intermediates,
                effective_n=args.effective_n,
                max_causal=args.max_causal,
                max_iter=args.max_iter,
            )
            analyzed.append(locus_report)

        aggregate_method_results(
            output_root=output_dir,
            scripts_dir=scripts_dir,
            input_root=rds_stage_dir,
            allow_empty=True,
        )

    n_analyzed = sum(1 for r in analyzed if r.get("status") == "analyzed")
    total_cs_rss = sum(
        int(r.get("methods", {}).get("susie_rss", {}).get("credible_sets") or 0)
        for r in analyzed
        if r.get("status") == "analyzed"
    )
    total_cs_rrss = sum(
        int(r.get("methods", {}).get("susie_rrss", {}).get("credible_sets") or 0)
        for r in analyzed
        if r.get("status") == "analyzed"
    )
    nonconverged = [
        r["region_id"]
        for r in analyzed
        if r.get("status") == "analyzed"
        and (
            not r.get("methods", {}).get("susie_rss", {}).get("converged", True)
            or not r.get("methods", {}).get("susie_rrss", {}).get("converged", True)
        )
    ]
    mismatch_flagged = [
        r["region_id"]
        for r in analyzed
        if r.get("status") == "analyzed"
        and (
            r.get("methods", {}).get("susie_rrss", {}).get("SER_fallback") is True
            or r.get("methods", {}).get("susie_rrss", {}).get("R_reliability_flag") is True
            or r.get("methods", {}).get("susie_rrss", {}).get("artifact_flag") is True
        )
    ]

    if n_analyzed == 0:
        overall_status = "NO_SIGNIFICANT_LOCI" if not candidate_regions else "NO_PANEL_MATCHED_LOCI"
    elif nonconverged or mismatch_flagged:
        overall_status = "WARN_DIAGNOSTIC_FLAGS"
    else:
        overall_status = "PASS"

    sample_sizes = [row["SS"] for row in significant if row["SS"] is not None]
    runtime_seconds = round(time.monotonic() - t_run0, 3)

    provenance = {
        "schema_version": 2,
        "study_id": study_id,
        "gwas_id": gwas_id,
        "trait": metadata.get("phenotype"),
        "gwas_ancestry": metadata.get("ancestry"),
        "ancestry_ld_status": ancestry_ld_status,
        "gwas_build": "GRCh38",
        "vcf_optional_tags": vcf_tags,
        "ld_reference": "Pan-UKB public AWS S3" if local_ld_panel is None else f"local_fixture:{local_ld_panel.name}",
        "ld_population": args.population,
        "ld_reference_sample_size": reference_n,
        "reference_sample_size_source": (
            f"UKBB.{args.population}.ldadj.variant.ht global n_samples"
            if args.reference_sample_size is None
            else "explicit_run_override"
        ),
        "liftover_chain": str(chain),
        "liftover_chain_sha256": sha256_file(chain),
        "target_reference_fasta": str(fasta),
        "target_reference_fai": str(fasta_index),
        "gene_resolution": gene_resolution,
        "gene_annotation": gene,
        "cis_window_bp": args.gene_window_bp,
        "significance_threshold_p": 10 ** (-args.significance_lp),
        "significance_threshold_lp": args.significance_lp,
        "cis_significance_threshold_lp": (
            args.cis_significance_lp
            if args.cis_significance_lp is not None
            else args.significance_lp
        ),
        "fine_mapping_window_definition": "kanai_greedy_lead_plus_minus_flank_merged",
        "region_flank_bp": args.region_flank_bp,
        "nominal_locus_window_bp": 2 * args.region_flank_bp,
        "max_region_width_bp": args.max_region_width_bp,
        "mhc_policy": "excluded_by_option" if args.exclude_mhc else "included_single_region",
        "mhc_region_grch38": f"{MHC_CHROM}:{MHC_START_GRCH38}-{MHC_END_GRCH38}",
        "mhc_core_grc_grch38": f"{MHC_CHROM}:{MHC_CORE_GRC_GRCH38[0]}-{MHC_CORE_GRC_GRCH38[1]}",
        "mhc_core_grc_grch37": f"{MHC_CHROM}:{MHC_CORE_GRC_GRCH37[0]}-{MHC_CORE_GRC_GRCH37[1]}",
        "max_loci": args.max_loci,
        "max_causal_L": args.max_causal,
        "max_iter": args.max_iter,
        "effective_n_override": args.effective_n,
        "keep_intermediates": args.keep_intermediates,
        "fine_mapping_methods": METHODS.split(","),
        "output_contract": "standalone_xz_rda_json_and_self_contained_html",
        "duplicate_panel_variant_policy": (
            "retain highest-LP GWAS row per Pan-UKB variant index; record exclusion counts in qc.json"
        ),
        "runtime_seconds": runtime_seconds,
    }

    finemapping_qc = {
        "schema_version": 2,
        "overall_status": overall_status,
        "study_id": study_id,
        "gwas_id": gwas_id,
        "ancestry_ld_status": ancestry_ld_status,
        "upstream_qc": qc_upstream,
        "genomewide_significant_variants": len(significant),
        "minimum_variant_ss_over_significant_variants": (
            min(sample_sizes) if sample_sizes else None
        ),
        "candidate_regions_total": len(candidate_regions),
        "selected_regions_count": len(selected_regions),
        "analyzed_regions_count": n_analyzed,
        "total_credible_sets_susie_rss": total_cs_rss,
        "total_credible_sets_susie_rrss": total_cs_rrss,
        "nonconverged_regions": nonconverged,
        "mismatch_flagged_regions": mismatch_flagged,
        "cis_gene": cis_summary,
        "mhc": mhc_summary,
        "loci": analyzed,
    }

    metadata_out = {
        **metadata,
        "finemapping_status": overall_status,
        "ld_population": args.population,
        "analyzed_regions_count": n_analyzed,
        "credible_sets_susie_rss": total_cs_rss,
        "credible_sets_susie_rrss": total_cs_rrss,
    }

    write_json(output_dir / "provenance.json", provenance)
    write_json(output_dir / "metadata.json", metadata_out)
    write_json(output_dir / "qc.json", finemapping_qc)

    if not args.no_report:
        render_gwas_report(output_dir, scripts_dir)

    if getattr(args, "output_record", None):
        write_json(Path(args.output_record), {
            "ok": True,
            "row_key": f"{study_id}/{gwas_id}",
            "study_id": study_id,
            "gwas_id": gwas_id,
            "settings_hash": getattr(args, "settings_hash", ""),
            "overall_status": overall_status,
            "analyzed_regions": n_analyzed,
            "cs_rss": total_cs_rss,
            "cs_rrss": total_cs_rrss,
            "output_dir": str(output_dir),
        })

    return output_dir


def _apply_row_json(args: argparse.Namespace) -> None:
    if not args.row_json:
        return
    row = json.loads(Path(args.row_json).read_text(encoding="utf-8"))
    args.vcf = row["vcf"]
    args.metadata = row.get("metadata") or None
    args.qc = row.get("qc") or None
    args.study_id = row["study_id"]
    args.gwas_id = row["gwas_id"]
    args.trait = row.get("trait") or None
    args.ancestry = row.get("ancestry") or None
    args.gene_id = row.get("gene_id") or None
    args.chain = row["chain"]
    args.grch37_fasta = row["grch37_fasta"]
    args.gene_annotation = row.get("gene_annotation") or None
    args.output_root = row["outdir"]
    args.cache_dir = row["cache_dir"]
    args.panukb_dir = row.get("panukb_dir") or None
    args.local_ld_panel = row.get("local_ld_panel") or None
    args.population = row["ld_population"]
    args.reference_sample_size = row.get("reference_sample_size")
    args.effective_n = row.get("effective_n")
    args.significance_lp = float(row["significance_lp"])
    args.cis_significance_lp = (
        float(row["cis_significance_lp"])
        if row.get("cis_significance_lp") is not None
        else None
    )
    args.region_flank_bp = int(row["region_flank_bp"])
    args.max_region_width_bp = int(row["max_region_width_bp"])
    args.gene_window_bp = int(row["gene_window_bp"])
    args.exclude_mhc = bool(row.get("exclude_mhc", False))
    args.max_loci = int(row.get("max_loci", 0))
    args.max_causal = int(row.get("max_causal", 10))
    args.max_iter = int(row.get("max_iter", 200))
    args.keep_intermediates = bool(row.get("keep_intermediates", False))
    args.force = bool(row.get("force", False))
    args.no_report = not bool(row.get("report", True))
    args.allow_ancestry_mismatch = bool(row.get("allow_ancestry_mismatch", False))
    args.settings_hash = row.get("settings_hash", "")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--row-json", help="JSON row specification emitted by prepare_manifest.py")
    parser.add_argument("--output-record", help="Path to write per-row completion/error JSON record")
    parser.add_argument("--vcf", help="Bgzipped + indexed GRCh38 GWAS-VCF (.vcf.gz)")
    parser.add_argument("--metadata", help="Optional GWASprocessing metadata.json")
    parser.add_argument("--qc", help="Optional GWASprocessing qc.json")
    parser.add_argument("--study-id", help="Study identifier (overrides or replaces metadata.json)")
    parser.add_argument("--gwas-id", help="GWAS identifier (overrides or replaces metadata.json)")
    parser.add_argument("--trait", help="Optional phenotype label")
    parser.add_argument("--ancestry", help="Optional GWAS ancestry label (e.g. EUR)")
    parser.add_argument("--gene-id", help="Optional explicit UniProt accession or HGNC gene symbol for cis_gene")
    parser.add_argument("--chain", help="GRCh38-to-GRCh37 UCSC chain")
    parser.add_argument("--grch37-fasta", help="GRCh37 primary assembly FASTA")
    parser.add_argument("--gene-annotation", help="GRCh38 UniProt/gene table (gene_regions_GRCh38.tsv)")
    parser.add_argument("--output-root", help="Root directory for <study_id>/<gwas_id>/ outputs")
    parser.add_argument("--cache-dir", default="~/.cache/gwasfinemapping")
    parser.add_argument(
        "--panukb-dir",
        help="Folder with installed whole-population indexes UKBB.<POP>.variant_index.parquet "
             "(default: panukb/ under GWASFM_REFERENCE_ROOT)",
    )
    parser.add_argument(
        "--local-ld-panel",
        help="Optional local .npz LD fixture for offline testing/benchmarking (replaces S3 ldcov call)",
    )
    parser.add_argument("--population", choices=tuple(REFERENCE_SAMPLES), default="EUR")
    parser.add_argument("--reference-sample-size", type=int)
    parser.add_argument("--effective-n", type=float, help="Optional explicit sample-size override for SuSiE")
    parser.add_argument("--significance-lp", type=float, default=-math.log10(5e-8))
    parser.add_argument(
        "--cis-significance-lp",
        type=float,
        help="Optional cis_gene significance threshold in -log10(P) (defaults to --significance-lp, P <= 5e-8)",
    )
    parser.add_argument(
        "--region-flank-bp",
        type=int,
        default=1_500_000,
        help="Half-width around greedy lead variants in bp (default: 1,500,000 = 3 Mb window, Kanai et al. 2021)",
    )
    parser.add_argument(
        "--max-region-width-bp",
        type=int,
        default=9_500_000,
        help="Maximum merged non-MHC region width in bp (<10 Mb Pan-UKB LD radius; default: 9,500,000)",
    )
    parser.add_argument(
        "--merge-distance-bp",
        type=int,
        default=250_000,
        help="Retained legacy option; greedy Kanai ±flank merging supersedes gap clustering.",
    )
    parser.add_argument("--gene-window-bp", type=int, default=1_000_000)
    parser.add_argument("--exclude-mhc", action="store_true", help="Exclude MHC (chr6:25.5-34.0 Mb GRCh38)")
    parser.add_argument("--max-loci", type=int, default=0, help="Maximum candidate regions to analyze; 0 = all")
    parser.add_argument("--max-causal", type=int, default=10, help="Maximum causal signals L in SuSiE (default: 10)")
    parser.add_argument("--max-iter", type=int, default=200, help="Maximum IBSS iterations in SuSiE (default: 200)")
    parser.add_argument("--keep-intermediates", action="store_true", help="Keep scratch LD/VCF/TSV files for debugging")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing GWAS output folder")
    parser.add_argument("--allow-ancestry-mismatch", action="store_true", help="Allow mismatched single-ancestry LD")
    parser.add_argument("--report", choices=("true", "false"), help="Whether to write per-GWAS report.html")
    parser.add_argument("--no-report", action="store_true", help="Skip HTML report rendering")

    args = parser.parse_args()
    args.settings_hash = ""
    if args.report is not None:
        args.no_report = args.report == "false"
    if args.row_json:
        _apply_row_json(args)
    if not all((args.vcf, args.chain, args.grch37_fasta, args.output_root)):
        parser.error("--vcf, --chain, --grch37-fasta, and --output-root are required (unless --row-json is supplied).")
    if args.max_loci < 0 or args.region_flank_bp < 1 or args.gene_window_bp < 0:
        parser.error("Locus count must be non-negative and window sizes valid.")
    if args.max_region_width_bp < 2 * args.region_flank_bp or args.max_region_width_bp > 9_900_000:
        parser.error("--max-region-width-bp must be >= 2 * --region-flank-bp and <= 9,900,000 (Pan-UKB 10 Mb limit).")
    if args.significance_lp <= 0 or (args.cis_significance_lp is not None and args.cis_significance_lp <= 0):
        parser.error("Significance thresholds (--significance-lp / --cis-significance-lp) must be positive.")

    if args.output_record:
        try:
            output = run(args)
            print(output)
        except Exception as exc:  # noqa: BLE001 - isolate row failure in batch mode & log error
            print(f"ERROR [{args.study_id}/{args.gwas_id}]: {exc}", file=sys.stderr)
            write_json(Path(args.output_record), {
                "ok": False,
                "row_key": f"{args.study_id}/{args.gwas_id}",
                "study_id": args.study_id or "",
                "gwas_id": args.gwas_id or "",
                "settings_hash": args.settings_hash,
                "error": str(exc),
            })
    else:
        output = run(args)
        print(output)


if __name__ == "__main__":
    main()
