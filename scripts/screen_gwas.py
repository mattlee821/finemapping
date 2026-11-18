#!/usr/bin/env python3
"""
Standalone Stage 1: Genome-wide significance screen (P <= 5e-8) & locus builder.

Can be run:
  1. Interactively / without arguments (uses the default local NPPB paths below), OR
  2. From CLI / Nextflow:
     python3 scripts/screen_gwas.py \
       --vcf /path/to/trait.gwas.vcf.gz \
       --metadata /path/to/trait.metadata.json \
       --qc /path/to/trait.qc.json \
       --out /path/to/screen_summary.json
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_finemapping import (  # noqa: E402
    REFERENCE_SAMPLES,
    build_analysis_regions,
    chromosome_sort_key,
    inspect_vcf_tags,
    load_or_init_sidecars,
    query_vcf,
    resolve_reference_sample_size,
    resolve_target_gene,
    write_json,
)

# ==============================================================================
# 1. Interactive defaults (used when run with no CLI flags)
# ==============================================================================
HOME = Path.home()
DEFAULT_NPPB_DIR = HOME / "research-local/GWASprocessing/data/GWAS/sun_2023_37794186/standardise/P16860_OID20049"
DEFAULT_REF_ROOT = Path(os.environ.get("GWASFM_REFERENCE_ROOT", str(HOME / ".cache/gwasfinemapping/reference")))
DEFAULT_OUT_JSON = HOME / "research-local/GWASfinemapping/test_nppb/screen_summary.json"


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
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", default=str(DEFAULT_NPPB_DIR / "P16860_OID20049.gwas.vcf.gz"))
    p.add_argument("--metadata", default=str(DEFAULT_NPPB_DIR / "P16860_OID20049.metadata.json"))
    p.add_argument("--qc", default=str(DEFAULT_NPPB_DIR / "P16860_OID20049.qc.json"))
    p.add_argument("--study-id", default="sun_2023_37794186")
    p.add_argument("--gwas-id", default="P16860_OID20049")
    p.add_argument("--trait", default=None)
    p.add_argument("--ancestry", default=None)
    p.add_argument("--gene-id", default=None)
    p.add_argument("--chain", default=str(DEFAULT_REF_ROOT / "genome/hg38ToHg19.over.chain.gz"))
    p.add_argument("--gene-annotation", default=str(DEFAULT_REF_ROOT / "annotation/gene_regions_GRCh38.tsv"))
    p.add_argument("--population", default="EUR")
    p.add_argument("--reference-sample-size", type=int, default=None)
    p.add_argument("--significance-lp", type=float, default=7.301029995664)
    p.add_argument("--cis-significance-lp", type=float, default=None)
    p.add_argument("--region-flank-bp", type=int, default=1_500_000)
    p.add_argument("--gene-window-bp", type=int, default=1_000_000)
    p.add_argument("--max-region-width-bp", type=int, default=9_500_000)
    p.add_argument("--exclude-mhc", action="store_true")
    p.add_argument("--allow-ancestry-mismatch", action="store_true")
    p.add_argument("--max-loci", type=int, default=2)
    p.add_argument("--out", default=str(DEFAULT_OUT_JSON))
    return p.parse_args(argv)


def screen_gwas(args: argparse.Namespace) -> dict[str, Any]:
    t0 = time.monotonic()
    u0, s0 = current_cpu_times()

    vcf = Path(args.vcf).expanduser().resolve()
    chain = Path(args.chain).expanduser().resolve()
    gene_path = Path(args.gene_annotation).expanduser().resolve() if args.gene_annotation else None

    if not vcf.is_file():
        raise FileNotFoundError(f"GWAS-VCF not found: {vcf}")
    if not (Path(f"{vcf}.tbi").is_file() or Path(f"{vcf}.csi").is_file()):
        raise FileNotFoundError(f"Tabix/CSI index missing for {vcf}")
    if not chain.is_file():
        raise FileNotFoundError(f"Chain file not found: {chain}")

    metadata, qc_upstream, study_id, gwas_id = load_or_init_sidecars(
        vcf=vcf,
        metadata_arg=args.metadata if args.metadata and Path(args.metadata).is_file() else None,
        qc_arg=args.qc if args.qc and Path(args.qc).is_file() else None,
        study_id_arg=args.study_id,
        gwas_id_arg=args.gwas_id,
        trait_arg=args.trait,
        ancestry_arg=args.ancestry,
    )
    gene, gene_resolution = resolve_target_gene(
        gene_path=gene_path if (gene_path and gene_path.is_file()) else None,
        gwas_id=gwas_id,
        metadata=metadata,
        explicit_gene_id=args.gene_id,
    )
    reference_n = resolve_reference_sample_size(args.population, args.reference_sample_size)
    if reference_n is None:
        raise ValueError(f"Supply --reference-sample-size for population {args.population}")

    gwas_ancestry = metadata.get("ancestry")
    if (
        gwas_ancestry in REFERENCE_SAMPLES
        and gwas_ancestry != args.population
        and not args.allow_ancestry_mismatch
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

    # --- Scan VCF for genome-wide significant variants ------------------------
    t_scan0 = time.monotonic()
    vcf_tags = inspect_vcf_tags(vcf)
    significant = query_vcf(vcf, minimum_lp=args.significance_lp, vcf_tags=vcf_tags)
    significant.sort(key=lambda r: (-r["LP"], chromosome_sort_key(r["CHROM"]), r["POS"]))
    scan_seconds = round(time.monotonic() - t_scan0, 3)

    # --- Build candidate regions (cis_gene, Kanai 3 Mb lead_locus, mhc) -------
    t_reg0 = time.monotonic()
    from pyliftover import LiftOver

    lifter = LiftOver(str(chain))
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
    region_build_seconds = round(time.monotonic() - t_reg0, 3)

    selected_regions = (
        candidate_regions[: args.max_loci] if args.max_loci > 0 else candidate_regions
    )
    min_sig_ss = min((r["SS"] for r in significant if r.get("SS") is not None), default=None)

    u1, s1 = current_cpu_times()
    wall_seconds = round(time.monotonic() - t0, 3)

    payload: dict[str, Any] = {
        "study_id": study_id,
        "gwas_id": gwas_id,
        "vcf": str(vcf),
        "vcf_tags": vcf_tags,
        "population": args.population,
        "reference_sample_size": reference_n,
        "ancestry_ld_status": ancestry_ld_status,
        "gene": gene,
        "gene_resolution": gene_resolution,
        "gene_window_bp": args.gene_window_bp,
        "significance_lp": args.significance_lp,
        "cis_significance_lp": args.cis_significance_lp,
        "region_flank_bp": args.region_flank_bp,
        "max_region_width_bp": args.max_region_width_bp,
        "exclude_mhc": args.exclude_mhc,
        "max_loci": args.max_loci,
        "genomewide_significant_variants": len(significant),
        "minimum_variant_ss_over_significant_variants": min_sig_ss,
        "status": "active" if len(selected_regions) > 0 else "no_active_regions",
        "cis_gene": cis_summary,
        "mhc": mhc_summary,
        "candidate_regions_total": len(candidate_regions),
        "selected_regions_count": len(selected_regions),
        "selected_regions": selected_regions,
        "metadata": metadata,
        "upstream_gwasprocessing_qc": qc_upstream,
        "timing": {
            "vcf_scan_seconds": scan_seconds,
            "region_build_seconds": region_build_seconds,
            "wall_seconds": wall_seconds,
            "user_seconds": round(u1 - u0, 3),
            "sys_seconds": round(s1 - s0, 3),
            "peak_rss_mb": current_peak_rss_mb(),
        },
    }

    out_path = Path(args.out).expanduser().resolve()
    write_json(out_path, payload)
    print(
        f"[screen_gwas] {study_id}/{gwas_id}: {len(significant)} sig variants -> "
        f"{len(candidate_regions)} candidate region(s), {len(selected_regions)} selected "
        f"({wall_seconds:.2f}s, peak RSS {payload['timing']['peak_rss_mb']:.1f} MB) -> {out_path}"
    )
    return payload


if __name__ == "__main__":
    screen_gwas(parse_args())
