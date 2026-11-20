#!/usr/bin/env python3
"""
Standalone Stage 4: Finalise one GWAS from screen_summary.json + per-locus prep & method outputs.

Combines per-locus region_result.rds files into GWAS-level method files (.rda and/or .rds),
writes {metadata,provenance,qc}.json (including per-step resource telemetry), renders
self-contained report.html, and cleans up scratch LD matrices unless --keep-intermediates is set.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import resource
import shutil
import sys
import time
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from run_finemapping import (  # noqa: E402
    MHC_CHROM,
    MHC_CORE_GRC_GRCH37,
    MHC_CORE_GRC_GRCH38,
    MHC_END_GRCH38,
    MHC_START_GRCH38,
    r_runtime_environment,
    run_checked,
    sha256_file,
    write_json,
)

HOME = Path.home()
DEFAULT_TEST_ROOT = HOME / "research-local/GWASfinemapping/test_nppb"
DEFAULT_REF_ROOT = Path(os.environ.get("GWASFM_REFERENCE_ROOT", str(HOME / ".cache/gwasfinemapping/reference")))


def normalize_method_id(name: str) -> str:
    val = name.strip().lower().replace("-", "_")
    if val == "polyfun_susie":
        return "polyfun"
    return val


def parse_methods_list(spec: str) -> list[str]:
    tokens = [normalize_method_id(t) for t in re.split(r"[,\s]+", spec.strip()) if t.strip()]
    seen: list[str] = []
    for t in tokens:
        if t not in seen:
            seen.append(t)
    return seen


def current_peak_rss_mb() -> float:
    ru_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ru_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak = max(ru_self, ru_child)
    if sys.platform == "darwin":
        return round(peak / (1024.0 * 1024.0), 1)
    return round(peak / 1024.0, 1)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--screen-json", default=str(DEFAULT_TEST_ROOT / "screen_summary.json"))
    p.add_argument("--loci-dir", default=str(DEFAULT_TEST_ROOT / "loci"))
    p.add_argument("--output-root", default=str(DEFAULT_TEST_ROOT / "published"))
    p.add_argument("--methods", default="p_ld,susie_rss,susie_rrss,carma,polyfun")
    p.add_argument("--output-format", default="rda", choices=["rda", "rds", "both"])
    p.add_argument("--compress", default="xz", choices=["xz", "gzip", "bzip2"])
    p.add_argument("--chain", default=str(DEFAULT_REF_ROOT / "genome/hg38ToHg19.over.chain.gz"))
    p.add_argument("--grch37-fasta", default=str(DEFAULT_REF_ROOT / "genome/GRCh37.fa"))
    p.add_argument("--max-causal", type=int, default=10)
    p.add_argument("--max-iter", type=int, default=100)
    p.add_argument("--cis-max-iter", type=int, default=200)
    p.add_argument("--min-maf", type=float, default=0.01)
    p.add_argument("--effective-n", type=float, default=None)
    p.add_argument("--report", default="true")
    p.add_argument("--keep-intermediates", action="store_true")
    p.add_argument("--cleanup-loci-scratch", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--settings-hash", default="")
    p.add_argument("--output-record", default=None)
    return p.parse_args(argv)


def finalise_gwas(args: argparse.Namespace) -> Path:
    t0 = time.monotonic()
    screen_path = Path(args.screen_json).expanduser().resolve()
    screen = json.loads(screen_path.read_text(encoding="utf-8"))

    study_id = screen["study_id"]
    gwas_id = screen["gwas_id"]
    metadata = dict(screen.get("metadata") or {})
    qc_upstream = screen.get("upstream_gwasprocessing_qc") or {}
    methods = parse_methods_list(args.methods)

    output_dir = Path(args.output_root).expanduser().resolve() / study_id / gwas_id
    loci_root = Path(args.loci_dir).expanduser().resolve() if args.loci_dir else None
    keep_names = {"logs"}
    if args.keep_intermediates or (loci_root is not None and loci_root.parent == output_dir):
        keep_names.add("loci")

    existing_non_logs = (
        [p for p in output_dir.iterdir() if p.name not in keep_names]
        if output_dir.exists()
        else []
    )
    if existing_non_logs:
        if args.force:
            for item in existing_non_logs:
                if item.is_dir():
                    shutil.rmtree(item)
                else:
                    item.unlink()
        else:
            raise FileExistsError(f"Output directory is not empty: {output_dir} (pass --force)")
    output_dir.mkdir(parents=True, exist_ok=True)

    selected_regions = screen.get("selected_regions") or []

    step_telemetry: list[dict[str, Any]] = []
    if "timing" in screen:
        step_telemetry.append({
            "step": "screen_gwas",
            "target": f"{study_id}/{gwas_id}",
            **screen["timing"],
        })

    analyzed_loci: list[dict[str, Any]] = []
    for reg in selected_regions:
        rid = reg["region_id"]
        loc_dir = loci_root / rid if loci_root else None
        prep_qc_file = loc_dir / "locus_prep_qc.json" if loc_dir else None
        if not prep_qc_file or not prep_qc_file.is_file():
            analyzed_loci.append({
                "region_id": rid,
                "region_type": reg["region_type"],
                "status": "missing_prep_output",
            })
            continue

        prep_qc = json.loads(prep_qc_file.read_text(encoding="utf-8"))
        if "timing" in prep_qc:
            step_telemetry.append({
                "step": "prep_locus",
                "target": rid,
                **prep_qc["timing"],
            })

        if prep_qc.get("status") != "prepared":
            analyzed_loci.append(prep_qc)
            continue

        method_qcs: dict[str, Any] = {}
        method_wall_sum = 0.0
        for m_id in methods:
            mqc_path = loc_dir / "methods" / m_id / "method_qc.json"
            if mqc_path.is_file():
                mqc = json.loads(mqc_path.read_text(encoding="utf-8"))
                t_m = mqc.pop("timing", {})
                method_qcs[m_id] = mqc
                method_wall_sum += float(t_m.get("wall_seconds", 0.0))
                step_telemetry.append({
                    "step": f"method_{m_id}",
                    "target": rid,
                    **t_m,
                })

        loc_entry = dict(prep_qc)
        loc_entry["status"] = "analyzed"
        loc_entry["methods"] = method_qcs
        loc_entry["wall_seconds"] = round(float(prep_qc.get("wall_seconds", 0.0)) + method_wall_sum, 3)
        analyzed_loci.append(loc_entry)

    failed_loci = [
        r for r in analyzed_loci
        if r.get("status") in ("missing_prep_output", "failed_locus")
        or (
            r.get("status") == "analyzed"
            and any(m not in (r.get("methods") or {}) for m in methods)
        )
    ]
    if failed_loci:
        # Remove only failed locus folders so already-completed sibling loci in `loci_root` are reused on retry!
        if not args.keep_intermediates and loci_root and loci_root.is_dir():
            for bad in failed_loci:
                bad_rid = bad.get("region_id")
                if bad_rid:
                    shutil.rmtree(loci_root / bad_rid, ignore_errors=True)
        bad_desc = ", ".join(
            f"{r['region_id']} ({r.get('error') or r.get('status')})" for r in failed_loci[:5]
        )
        raise RuntimeError(f"GWAS {study_id}/{gwas_id} had {len(failed_loci)} failed locus/loci: {bad_desc}")

    # --- Aggregate per-region .rds files into GWAS-level method files ---------
    t_agg0 = time.monotonic()
    r_env = r_runtime_environment(PROJECT_ROOT)
    agg_cmd = [
        "Rscript", str(SCRIPT_DIR / "aggregate_method_results.R"),
        "--output-root", str(output_dir),
        "--methods", ",".join(methods),
        "--format", args.output_format,
        "--compress", args.compress,
        "--allow-empty", "true",
        "--keep-rds", "true" if not args.cleanup_loci_scratch else "false",
    ]
    if loci_root and loci_root.is_dir():
        agg_cmd.extend(["--input-root", str(loci_root)])
    run_checked(agg_cmd, env=r_env, cwd=PROJECT_ROOT)
    agg_seconds = round(time.monotonic() - t_agg0, 3)
    step_telemetry.append({
        "step": "aggregate_methods",
        "target": f"{study_id}/{gwas_id}",
        "wall_seconds": agg_seconds,
    })

    # Clean up per-GWAS locus scratch immediately once all loci have been aggregated
    if args.cleanup_loci_scratch and not args.keep_intermediates and loci_root and loci_root.is_dir():
        shutil.rmtree(loci_root, ignore_errors=True)
        if loci_root.parent.parent.name == "by_gwas":
            shutil.rmtree(loci_root.parent, ignore_errors=True)

    n_analyzed = sum(1 for r in analyzed_loci if r.get("status") == "analyzed")
    total_cs_rss = sum(
        int(r.get("methods", {}).get("susie_rss", {}).get("credible_sets") or 0)
        for r in analyzed_loci if r.get("status") == "analyzed"
    )
    total_cs_rrss = sum(
        int(r.get("methods", {}).get("susie_rrss", {}).get("credible_sets") or 0)
        for r in analyzed_loci if r.get("status") == "analyzed"
    )
    total_cs_carma = sum(
        int(r.get("methods", {}).get("carma", {}).get("credible_sets") or 0)
        for r in analyzed_loci if r.get("status") == "analyzed"
    )
    total_cs_polyfun = sum(
        int(r.get("methods", {}).get("polyfun", {}).get("credible_sets") or 0)
        for r in analyzed_loci if r.get("status") == "analyzed"
    )

    nonconverged = [
        r["region_id"]
        for r in analyzed_loci
        if r.get("status") == "analyzed"
        and any(
            r.get("methods", {}).get(m, {}).get("converged") is False
            for m in ("susie_rss", "susie_rrss", "carma", "polyfun")
            if m in r.get("methods", {})
        )
    ]
    suppressed_cs_regions = [
        r["region_id"]
        for r in analyzed_loci
        if r.get("status") == "analyzed"
        and any(
            r.get("methods", {}).get(m, {}).get("credible_sets_suppressed") is True
            for m in ("susie_rss", "susie_rrss", "polyfun")
            if m in r.get("methods", {})
        )
    ]
    mismatch_flagged = [
        r["region_id"]
        for r in analyzed_loci
        if r.get("status") == "analyzed"
        and (
            r.get("methods", {}).get("susie_rrss", {}).get("SER_fallback") is True
            or r.get("methods", {}).get("susie_rrss", {}).get("R_reliability_flag") is True
            or r.get("methods", {}).get("susie_rrss", {}).get("artifact_flag") is True
            or int(r.get("methods", {}).get("carma", {}).get("outlier_count") or 0) > 0
        )
    ]

    if n_analyzed == 0:
        overall_status = (
            "NO_SIGNIFICANT_LOCI"
            if int(screen.get("candidate_regions_total", 0)) == 0
            else "NO_PANEL_MATCHED_LOCI"
        )
    elif nonconverged or suppressed_cs_regions or mismatch_flagged:
        overall_status = "WARN_DIAGNOSTIC_FLAGS"
    else:
        overall_status = "PASS"

    chain_p = Path(args.chain).expanduser().resolve()
    fasta_p = Path(args.grch37_fasta).expanduser().resolve()
    sig_lp = float(screen["significance_lp"])
    cis_sig_lp = screen.get("cis_significance_lp")
    flank_bp = int(screen["region_flank_bp"])

    provenance = {
        "schema_version": 3,
        "study_id": study_id,
        "gwas_id": gwas_id,
        "trait": metadata.get("phenotype"),
        "gwas_ancestry": metadata.get("ancestry"),
        "ancestry_ld_status": screen["ancestry_ld_status"],
        "gwas_build": "GRCh38",
        "vcf_optional_tags": screen.get("vcf_tags", {}),
        "ld_reference": "Pan-UKB public AWS S3",
        "ld_population": screen["population"],
        "ld_reference_sample_size": screen["reference_sample_size"],
        "liftover_chain": str(chain_p),
        "liftover_chain_sha256": sha256_file(chain_p) if chain_p.is_file() else "",
        "target_reference_fasta": str(fasta_p),
        "target_reference_fai": f"{fasta_p}.fai",
        "gene_resolution": screen.get("gene_resolution"),
        "gene_annotation": screen.get("gene"),
        "cis_window_bp": screen.get("gene_window_bp", 1_000_000),
        "significance_threshold_p": 10 ** (-sig_lp),
        "significance_threshold_lp": sig_lp,
        "cis_significance_threshold_lp": cis_sig_lp if cis_sig_lp is not None else sig_lp,
        "min_maf": float(screen.get("min_maf", args.min_maf)),
        "fine_mapping_window_definition": "kanai_greedy_lead_plus_minus_flank_merged",
        "region_flank_bp": flank_bp,
        "nominal_locus_window_bp": 2 * flank_bp,
        "max_region_width_bp": screen.get("max_region_width_bp", 9_500_000),
        "mhc_policy": "excluded_by_option" if screen.get("exclude_mhc") else "included_single_region",
        "mhc_region_grch38": f"{MHC_CHROM}:{MHC_START_GRCH38}-{MHC_END_GRCH38}",
        "mhc_core_grc_grch38": f"{MHC_CHROM}:{MHC_CORE_GRC_GRCH38[0]}-{MHC_CORE_GRC_GRCH38[1]}",
        "mhc_core_grc_grch37": f"{MHC_CHROM}:{MHC_CORE_GRC_GRCH37[0]}-{MHC_CORE_GRC_GRCH37[1]}",
        "max_loci": screen.get("max_loci", 0),
        "max_causal_L": args.max_causal,
        "max_iter": args.max_iter,
        "cis_max_iter": args.cis_max_iter,
        "effective_n_override": args.effective_n,
        "keep_intermediates": args.keep_intermediates,
        "fine_mapping_methods": methods,
        "output_format": args.output_format,
        "output_compression": args.compress,
    }

    finemapping_qc = {
        "schema_version": 3,
        "overall_status": overall_status,
        "study_id": study_id,
        "gwas_id": gwas_id,
        "ancestry_ld_status": screen["ancestry_ld_status"],
        "upstream_qc": qc_upstream,
        "genomewide_significant_variants": screen["genomewide_significant_variants"],
        "significant_variants_removed_min_maf": screen.get("significant_variants_removed_min_maf", 0),
        "total_variants_removed_min_maf": sum(int(r.get("variants_removed_min_maf") or 0) for r in analyzed_loci),
        "minimum_variant_ss_over_significant_variants": screen.get(
            "minimum_variant_ss_over_significant_variants"
        ),
        "candidate_regions_total": screen["candidate_regions_total"],
        "selected_regions_count": screen["selected_regions_count"],
        "analyzed_regions_count": n_analyzed,
        "total_credible_sets_susie_rss": total_cs_rss,
        "total_credible_sets_susie_rrss": total_cs_rrss,
        "total_credible_sets_carma": total_cs_carma,
        "total_credible_sets_polyfun": total_cs_polyfun,
        "nonconverged_regions": nonconverged,
        "suppressed_credible_set_regions": suppressed_cs_regions,
        "mismatch_flagged_regions": mismatch_flagged,
        "cis_gene": screen.get("cis_gene"),
        "mhc": screen.get("mhc"),
        "loci": analyzed_loci,
        "step_telemetry": step_telemetry,
    }

    metadata_out = {
        **metadata,
        "finemapping_status": overall_status,
        "ld_population": screen["population"],
        "fine_mapping_methods": methods,
        "analyzed_regions_count": n_analyzed,
        "credible_sets_susie_rss": total_cs_rss,
        "credible_sets_susie_rrss": total_cs_rrss,
        "credible_sets_carma": total_cs_carma,
        "credible_sets_polyfun": total_cs_polyfun,
    }

    write_json(output_dir / "provenance.json", provenance)
    write_json(output_dir / "metadata.json", metadata_out)
    write_json(output_dir / "qc.json", finemapping_qc)

    do_report = str(args.report).lower() in ("1", "true", "yes", "y")
    if do_report:
        t_rep0 = time.monotonic()
        run_checked(
            ["Rscript", str(SCRIPT_DIR / "render_report.R"), "--gwas-dir", str(output_dir)],
            env=r_env,
            cwd=PROJECT_ROOT,
        )
        rep_seconds = round(time.monotonic() - t_rep0, 3)
        step_telemetry.append({
            "step": "render_report",
            "target": f"{study_id}/{gwas_id}",
            "wall_seconds": rep_seconds,
        })
        finemapping_qc["step_telemetry"] = step_telemetry
        write_json(output_dir / "qc.json", finemapping_qc)

    provenance["runtime_seconds"] = round(time.monotonic() - t0, 3)
    provenance["finalise_peak_rss_mb"] = current_peak_rss_mb()
    write_json(output_dir / "provenance.json", provenance)

    if args.output_record:
        write_json(Path(args.output_record).expanduser().resolve(), {
            "ok": True,
            "row_key": f"{study_id}/{gwas_id}",
            "study_id": study_id,
            "gwas_id": gwas_id,
            "settings_hash": args.settings_hash,
            "overall_status": overall_status,
            "analyzed_regions": n_analyzed,
            "cs_rss": total_cs_rss,
            "cs_rrss": total_cs_rrss,
            "output_dir": str(output_dir),
        })

    print(
        f"[finalise_gwas] {study_id}/{gwas_id}: status={overall_status} | "
        f"regions={n_analyzed} | methods={','.join(methods)} -> {output_dir}"
    )
    return output_dir


if __name__ == "__main__":
    finalise_gwas(parse_args())
