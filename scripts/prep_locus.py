#!/usr/bin/env python3
"""
Standalone Stage 2: Shared per-locus data preparation (liftover, panel match, LD extraction).

Produces shared inputs consumed by all standalone fine-mapping method scripts:
  - <locus_dir>/matched_variants.tsv   (aligned GRCh38 + GRCh37 table with EZ and LD_R2_TO_LEAD)
  - <locus_dir>/matrix.npy             (Fortran-order float64 normalized symmetric LD matrix)
  - <locus_dir>/matrix.ids.txt         (1-line-per-variant companion IDs)
  - <locus_dir>/locus_prep_qc.json     (variant attrition, LD validation, SS stats, sub-step timings)

Can be run:
  1. Interactively / with no arguments (uses the default local NPPB cis_NPPB locus & local test BM if present), OR
  2. From CLI / Nextflow:
     python3 scripts/prep_locus.py \
       --screen-json /path/to/screen_summary.json \
       --locus-index 0 \
       --outdir /path/to/loci/cis_NPPB_chr1_10857464_12858945
"""

from __future__ import annotations

import argparse
from collections import Counter
import gc
import json
import os
from pathlib import Path
import resource
import sys
import tempfile
import time
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_finemapping import (  # noqa: E402
    MHC_CHROM,
    MHC_CORE_GRC_GRCH37,
    _normalize_ld_matrix_inplace,
    _write_symmetric_fortran_npy,
    chromosome_sort_key,
    default_panukb_dir,
    ensure_regional_index,
    extract_ld_from_local_panel,
    format_vcf_region,
    inspect_vcf_tags,
    is_mhc_core_grc_grch38,
    is_mhc_grch38,
    match_panel_variants,
    matrix_to_npy,
    normalize_chromosome,
    normalize_lifted_vcf,
    query_vcf,
    resolve_panel_matches,
    write_json,
    write_lifted_vcf,
    write_tsv,
)

# ==============================================================================
# 1. Interactive defaults (used when run without CLI arguments)
# ==============================================================================
HOME = Path.home()
DEFAULT_SCREEN_JSON = HOME / "research-local/GWASfinemapping/test_nppb/screen_summary.json"
DEFAULT_REF_ROOT = Path(os.environ.get("GWASFM_REFERENCE_ROOT", str(HOME / ".cache/gwasfinemapping/reference")))
DEFAULT_CACHE_DIR = Path(os.environ.get("GWASFM_CACHE_DIR", str(HOME / ".cache/gwasfinemapping")))
DEFAULT_LOCAL_TEST_BM = HOME / "research-local/GWASfinemapping/nppb_test_bm/UKBB.EUR.ldadj.bm"
DEFAULT_LOCUS_ROOT = HOME / "research-local/GWASfinemapping/test_nppb/loci"


def current_peak_rss_mb() -> float:
    ru_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ru_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak = max(ru_self, ru_child)
    if sys.platform == "darwin":
        return round(peak / (1024.0 * 1024.0), 1)
    return round(peak / 1024.0, 1)


def current_cpu_times() -> tuple[float, float]:
    s = resource.getrusage(resource.RUSAGE_SELF)
    c = resource.getrusage(resource.RUSAGE_CHILDREN)
    return round(s.ru_utime + c.ru_utime, 3), round(s.ru_stime + c.ru_stime, 3)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    no_cli = (argv is None and len(sys.argv) == 1)
    default_bm = (
        str(DEFAULT_LOCAL_TEST_BM)
        if (no_cli and DEFAULT_LOCAL_TEST_BM.is_dir())
        else None
    )
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--screen-json", default=str(DEFAULT_SCREEN_JSON))
    p.add_argument("--locus-index", type=int, default=0, help="0-based index in screen_summary.selected_regions")
    p.add_argument("--region-id", default=None, help="Select locus by region_id instead of --locus-index")
    p.add_argument("--outdir", default=None, help="Output locus folder (defaults to <test_nppb>/loci/<region_id>)")
    p.add_argument("--chain", default=str(DEFAULT_REF_ROOT / "genome/hg38ToHg19.over.chain.gz"))
    p.add_argument("--grch37-fasta", default=str(DEFAULT_REF_ROOT / "genome/GRCh37.fa"))
    p.add_argument("--panukb-dir", default=None)
    p.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR))
    p.add_argument("--bm-path", default=default_bm, help="Optional BlockMatrix path override (defaults to live S3)")
    p.add_argument("--local-ld-panel", default=None, help="Optional offline .npz fixture LD panel")
    p.add_argument("--fetch-workers", type=int, default=8)
    p.add_argument("--block-cache", type=int, default=4)
    p.add_argument("--keep-intermediates", action="store_true")
    return p.parse_args(argv)


def prep_locus(args: argparse.Namespace) -> dict[str, Any]:
    t0 = time.monotonic()
    u0, s0 = current_cpu_times()

    screen_path = Path(args.screen_json).expanduser().resolve()
    if not screen_path.is_file():
        raise FileNotFoundError(
            f"screen_summary.json not found at {screen_path} (run scripts/screen_gwas.py first)."
        )
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    regions = screen.get("selected_regions") or []
    if not regions:
        raise ValueError(f"No active regions in {screen_path}")

    if args.region_id:
        matches_reg = [r for r in regions if r["region_id"] == args.region_id]
        if not matches_reg:
            raise ValueError(f"Region {args.region_id} not found in {screen_path}")
        locus = matches_reg[0]
    else:
        if args.locus_index < 0 or args.locus_index >= len(regions):
            raise IndexError(f"locus_index {args.locus_index} out of range (0..{len(regions)-1})")
        locus = regions[args.locus_index]

    region_id = locus["region_id"]
    region_type = locus["region_type"]
    chrom = normalize_chromosome(locus["CHROM"])
    start = int(locus["start"])
    end = int(locus["end"])

    outdir = (
        Path(args.outdir).expanduser().resolve()
        if args.outdir
        else (DEFAULT_LOCUS_ROOT / region_id).resolve()
    )
    outdir.mkdir(parents=True, exist_ok=True)

    vcf = Path(screen["vcf"]).resolve()
    vcf_tags = screen.get("vcf_tags") or inspect_vcf_tags(vcf)
    population = screen["population"]
    gene = screen.get("gene")
    gene_window_bp = int(screen.get("gene_window_bp", 1_000_000))
    chain = Path(args.chain).expanduser().resolve()
    fasta = Path(args.grch37_fasta).expanduser().resolve()
    cache_dir = Path(args.cache_dir).expanduser().resolve()
    panukb_dir = (
        Path(args.panukb_dir).expanduser().resolve()
        if args.panukb_dir
        else default_panukb_dir()
    )
    local_ld_panel = Path(args.local_ld_panel).expanduser().resolve() if args.local_ld_panel else None

    # --- 2. Regional VCF extraction -------------------------------------------
    t_vcf0 = time.monotonic()
    records = query_vcf(
        vcf,
        region=format_vcf_region(vcf_tags, chrom, start, end),
        vcf_tags=vcf_tags,
    )
    vcf_query_seconds = round(time.monotonic() - t_vcf0, 3)

    # --- 3. GRCh38 -> GRCh37 liftover + FASTA-verified normalisation ----------
    t_lift0 = time.monotonic()
    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
    with tempfile.TemporaryDirectory(prefix=f"gwasfm-lift-{region_id[:20]}-") as lift_tmp:
        lifted_vcf_path = outdir / "lifted.grch37.vcf" if args.keep_intermediates else Path(lift_tmp) / "lifted.vcf"
        lifted_by_id, liftover_excluded = write_lifted_vcf(records, lifter, lifted_vcf_path)
        norm_tsv_out = outdir / "normalized.grch37.tsv" if args.keep_intermediates else None
        normalized = normalize_lifted_vcf(lifted_vcf_path, fasta, norm_tsv_out) if lifted_by_id else []
    liftover_seconds = round(time.monotonic() - t_lift0, 3)

    for row in normalized:
        row.update({
            key: value
            for key, value in lifted_by_id[row["source_variant_key"]].items()
            if key not in ("LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "ES", "SE", "LP", "SS")
        })

    cross_chrom_excluded = [r for r in normalized if r["LD_CHROM"] != chrom]
    for r in cross_chrom_excluded:
        r["status"] = "excluded"
        r["reason"] = "lifted_to_different_chromosome"
    normalized = [r for r in normalized if r["LD_CHROM"] == chrom]

    if not normalized:
        qc_skip = {
            "region_id": region_id,
            "region_type": region_type,
            "source_region_grch38": f"{chrom}:{start}-{end}",
            "status": "skipped_no_lifted_variants",
            "mhc_region": locus["mhc_region"],
            "mhc_core_grc": locus["mhc_core_grc"],
            "variants_in_source_region": len(records),
            "variants_lifted": 0,
            "variants_ld_matched": 0,
            "exclusion_reasons": dict(Counter(r["reason"] for r in liftover_excluded + cross_chrom_excluded)),
        }
        write_json(outdir / "locus_prep_qc.json", qc_skip)
        return qc_skip

    target_chrom = chrom
    target_start = max(1, min(r["LD_POS"] for r in normalized) - 100)
    target_end = max(r["LD_POS"] for r in normalized) + 100
    grch37_span = target_end - target_start
    if grch37_span > 9_990_000:
        raise ValueError(
            f"Region {region_id} spans {grch37_span:,} bp on GRCh37, exceeding Pan-UKB's 10 Mb LD band."
        )

    # --- 4. Regional Pan-UKB VariantIndex slice & vectorized allele match -----
    t_idx0 = time.monotonic()
    index_path = ensure_regional_index(
        population, target_chrom, target_start, target_end, cache_dir, SCRIPT_DIR, panukb_dir
    )
    index_cut_seconds = round(time.monotonic() - t_idx0, 3)

    t_m0 = time.monotonic()
    matches = match_panel_variants(index_path, normalized)
    matched, panel_excluded = resolve_panel_matches(normalized, matches)
    match_seconds = round(time.monotonic() - t_m0, 3)

    all_excluded = liftover_excluded + cross_chrom_excluded + panel_excluded
    reason_counts = dict(Counter(r["reason"] for r in all_excluded))

    if not matched:
        qc_skip = {
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
        write_json(outdir / "locus_prep_qc.json", qc_skip)
        return qc_skip

    matched.sort(key=lambda r: (chromosome_sort_key(r["LD_CHROM"]), r["LD_POS"], r["LD_REF"], r["LD_ALT"]))
    original_by_key = {(r["CHROM"], r["POS"], r["REF"], r["ALT"]): r for r in records}

    matched_rows: list[dict[str, Any]] = []
    for row in matched:
        source = original_by_key.get((row["CHROM"], row["POS"], row["REF"], row["ALT"]))
        if source is None:
            raise ValueError("Could not recover source GWAS variant row.")
        cis_trans_label = (
            "cis"
            if gene is not None
            and source["CHROM"] == gene["chrom"]
            and (gene["start"] - gene_window_bp) <= source["POS"] <= (gene["end"] + gene_window_bp)
            else ("trans" if gene is not None else "na")
        )
        gwas_variant_id = f"{source['CHROM']}:{source['POS']}:{source['REF']}:{source['ALT']}"
        matched_rows.append({
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
        })

    # --- 5. BlockMatrix LD extraction + in-place Fortran .npy write -----------
    t_ld0 = time.monotonic()
    matrix_npy_path = outdir / "matrix.npy"
    lead_idx = max(range(len(matched_rows)), key=lambda idx: (matched_rows[idx]["LP"], -matched_rows[idx]["GWAS_POS"]))

    if local_ld_panel is not None:
        import numpy as np

        ld_prefix = outdir / f"panukb_{population.lower()}"
        extract_ld_from_local_panel(local_ld_panel, matched_rows, ld_prefix)
        ld_validation = matrix_to_npy(
            Path(f"{ld_prefix}.npz"),
            Path(f"{ld_prefix}.variants.tsv"),
            matrix_npy_path,
        )
        arr = np.load(matrix_npy_path, mmap_mode="r")
        lead_r2 = (arr[:, lead_idx] ** 2).tolist()
        del arr
        if not args.keep_intermediates:
            Path(f"{ld_prefix}.npz").unlink(missing_ok=True)
            Path(f"{ld_prefix}.variants.tsv").unlink(missing_ok=True)
    else:
        from ldcov.io.blockmatrix.reader import HailBlockMatrixReader
        from ldcov.ld_bm.extract import _apply_flips, _assemble

        bm_uri = (
            args.bm_path
            if args.bm_path
            else f"s3://pan-ukb-us-east-1/ld_release/UKBB.{population}.ldadj.bm"
        )
        global_idxs = [int(r["PANEL_IDX"]) for r in matched_rows]
        reader = HailBlockMatrixReader(bm_uri, block_cache=args.block_cache)
        matrix = _assemble(reader, global_idxs, max_workers=args.fetch_workers)
        del reader
        flip_positions = [i for i, r in enumerate(matched_rows) if r["PANEL_ALLELE_SWAPPED"]]
        matrix = _apply_flips(matrix, flip_positions)

        ids = [
            f"{normalize_chromosome(str(r['LD_CHROM']))}:{int(r['LD_POS'])}:{r['LD_REF']}:{r['LD_ALT']}"
            for r in matched_rows
        ]
        ld_validation = _normalize_ld_matrix_inplace(matrix, ids)
        lead_r2 = (matrix[:, lead_idx] ** 2).tolist()
        _write_symmetric_fortran_npy(matrix, ids, matrix_npy_path)
        del matrix
        gc.collect()
    ld_extract_seconds = round(time.monotonic() - t_ld0, 3)

    for idx, r2_val in enumerate(lead_r2):
        matched_rows[idx]["LD_R2_TO_LEAD"] = r2_val

    matched_variants_path = outdir / "matched_variants.tsv"
    write_tsv(
        matched_variants_path,
        matched_rows,
        [
            "CHROM", "POS", "ID", "REF", "ALT",
            "GWAS_VARIANT_ID", "GWAS_CHROM", "GWAS_POS", "GWAS_REF", "GWAS_ALT",
            "LP", "ES", "SE", "SS", "AF", "AFS", "HM", "EZ", "LD_R2_TO_LEAD",
            "LD_CHROM", "LD_POS", "LD_REF", "LD_ALT", "LD_VARIANT_ID",
            "PANEL_IDX", "PANEL_ALLELE_SWAPPED", "cis_trans",
            "mhc_region", "mhc_core_grc",
        ],
    )

    n_swapped = sum(1 for r in matched_rows if r["PANEL_ALLELE_SWAPPED"])
    lead_matched = matched_rows[lead_idx]
    ss_vals = sorted(float(r["SS"]) for r in matched_rows if r.get("SS") is not None)
    ss_min = ss_vals[0] if ss_vals else None
    ss_max = ss_vals[-1] if ss_vals else None
    ss_median = ss_vals[len(ss_vals) // 2] if ss_vals else None
    ss_ratio = round(ss_min / ss_max, 4) if (ss_min is not None and ss_max and ss_max > 0) else None
    ss_het_warn = bool(ss_ratio is not None and ss_ratio < 0.80)

    u1, s1 = current_cpu_times()
    wall_seconds = round(time.monotonic() - t0, 3)
    npy_mb = round(matrix_npy_path.stat().st_size / (1024.0 * 1024.0), 2)

    prep_qc: dict[str, Any] = {
        "region_id": region_id,
        "region_type": region_type,
        "status": "prepared",
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
        "index_cut_seconds": index_cut_seconds,
        "wall_seconds": wall_seconds,
        "ld_validation": ld_validation,
        "bm_source": "local_override" if args.bm_path else ("local_npz" if local_ld_panel else "live_s3"),
        "scratch_matrix_mb": npy_mb,
        "timing": {
            "vcf_query_seconds": vcf_query_seconds,
            "liftover_norm_seconds": liftover_seconds,
            "index_cut_seconds": index_cut_seconds,
            "variant_match_seconds": match_seconds,
            "ld_extract_normalize_seconds": ld_extract_seconds,
            "wall_seconds": wall_seconds,
            "user_seconds": round(u1 - u0, 3),
            "sys_seconds": round(s1 - s0, 3),
            "peak_rss_mb": current_peak_rss_mb(),
        },
    }
    write_json(outdir / "locus_prep_qc.json", prep_qc)
    print(
        f"[prep_locus] {region_id}: {len(matched_rows)}/{len(records)} matched | "
        f"matrix={npy_mb:.1f} MB | ld_extract={ld_extract_seconds:.2f}s, total={wall_seconds:.2f}s, "
        f"peak RSS {prep_qc['timing']['peak_rss_mb']:.1f} MB -> {outdir}"
    )
    return prep_qc


if __name__ == "__main__":
    prep_locus(parse_args())
