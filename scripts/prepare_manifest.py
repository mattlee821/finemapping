#!/usr/bin/env python3
"""Prepare, preflight, state-track, and copy logs for GWASfinemapping runs
(mirrors GWASprocessing's prepare_manifest.py, update_state.py, and copy_logs.py)."""

from __future__ import annotations

import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

REFERENCE_POPULATIONS = {"AFR", "AMR", "CSA", "EAS", "EUR", "MID"}
METHOD_STEMS = {
    "p_ld": "p-ld",
    "susie_rss": "susie-rss",
    "susie_rrss": "susie-rrss",
    "carma": "carma",
    "polyfun": "polyfun",
}


def normalize_methods_spec(spec: str) -> list[str]:
    tokens = [
        t.strip().lower().replace("-", "_")
        for t in re.split(r"[,\s]+", (spec or "").strip())
        if t.strip()
    ]
    out: list[str] = []
    for t in tokens:
        if t == "polyfun_susie":
            t = "polyfun"
        if t not in METHOD_STEMS:
            raise ValueError(f"Unsupported fine-mapping method '{t}' (allowed: {sorted(METHOD_STEMS)})")
        if t not in out:
            out.append(t)
    if not out:
        raise ValueError("At least one fine-mapping method must be specified.")
    return out
STATE_COLUMNS = [
    "row_key",
    "study_id",
    "gwas_id",
    "settings_hash",
    "status",
    "analyzed_regions",
    "cs_rss",
    "cs_rrss",
    "output_dir",
]
RUN_MANIFEST_COLUMNS = [
    "row_id",
    "row_key",
    "study_id",
    "gwas_id",
    "trait",
    "ancestry",
    "gene_id",
    "ld_population",
    "vcf",
    "vcf_index",
    "metadata",
    "qc",
    "settings_hash",
]


def is_ignored_entry(path: Path) -> bool:
    return path.name.startswith(".") or path.name.startswith("._")


def normalized_build(value: Any) -> str:
    if value is None:
        return ""
    norm = re.sub(r"[^a-z0-9]", "", str(value).lower())
    if norm.startswith(("grch38", "hg38")):
        return "GRCh38"
    if norm.startswith(("grch37", "hg19")):
        return "GRCh37"
    return ""


def find_index_for_vcf(vcf: Path, explicit_index: str | None = None) -> Path | None:
    if explicit_index:
        p = Path(explicit_index).expanduser().resolve()
        return p if p.is_file() else None
    for ext in (".tbi", ".csi"):
        cand = Path(f"{vcf}{ext}")
        if cand.is_file():
            return cand
    return None


def resolve_gwas_folder(folder: Path, fallback_study_id: str | None = None) -> dict[str, Any] | None:
    """Inspect one directory to see if it holds a standardised GWAS-VCF (supports both
    GWASprocessing's `standardised.gwas.vcf.gz` layout and `<output_id>.gwas.vcf.gz`)."""
    if not folder.is_dir() or is_ignored_entry(folder):
        return None

    dirname = folder.name
    vcf_candidates = [
        folder / "standardised.gwas.vcf.gz",
        folder / f"{dirname}.gwas.vcf.gz",
        folder / f"{dirname}.vcf.gz",
    ]
    vcf = next((p for p in vcf_candidates if p.is_file()), None)
    if vcf is None:
        vcfs = [
            p for p in sorted(folder.glob("*.vcf.gz"))
            if not is_ignored_entry(p)
        ]
        if len(vcfs) == 1:
            vcf = vcfs[0]
    if vcf is None:
        return None

    meta_candidates = [
        folder / "standardised.metadata.json",
        folder / f"{dirname}.metadata.json",
        folder / "metadata.json",
    ]
    qc_candidates = [
        folder / "standardised.qc.json",
        folder / f"{dirname}.qc.json",
        folder / "qc.json",
    ]
    meta_path = next((p for p in meta_candidates if p.is_file()), None)
    qc_path = next((p for p in qc_candidates if p.is_file()), None)

    inferred_study = fallback_study_id
    if not inferred_study:
        if folder.parent.name == "standardise" and folder.parent.parent.name:
            inferred_study = folder.parent.parent.name
        else:
            inferred_study = folder.parent.name or "study"

    return {
        "study_id": inferred_study,
        "gwas_id": dirname,
        "vcf": str(vcf.resolve()),
        "vcf_index": str(find_index_for_vcf(vcf) or ""),
        "metadata": str(meta_path.resolve()) if meta_path else "",
        "qc": str(qc_path.resolve()) if qc_path else "",
        "trait": "",
        "ancestry": "",
        "gene_id": "",
    }


def discover_input_dir(input_dir: Path) -> list[dict[str, Any]]:
    root = input_dir.expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f"Input directory does not exist: {root}")

    direct = resolve_gwas_folder(root)
    if direct is not None:
        return [direct]

    discovered: list[dict[str, Any]] = []
    for child in sorted(root.iterdir(), key=lambda p: p.name):
        if is_ignored_entry(child) or not child.is_dir():
            continue
        as_gwas = resolve_gwas_folder(
            child,
            fallback_study_id=root.parent.name if root.name == "standardise" else root.name,
        )
        if as_gwas is not None:
            discovered.append(as_gwas)
            continue
        std_dir = child / "standardise"
        if std_dir.is_dir():
            for sub in sorted(std_dir.iterdir(), key=lambda p: p.name):
                if is_ignored_entry(sub) or not sub.is_dir():
                    continue
                sub_gwas = resolve_gwas_folder(sub, fallback_study_id=child.name)
                if sub_gwas is not None:
                    discovered.append(sub_gwas)
    return discovered


def detect_sidecars_near_vcf(vcf: Path) -> tuple[str, str]:
    folder = vcf.parent
    stem = re.sub(r"\.(gwas\.)?vcf(\.gz|\.bgz)?$", "", vcf.name)
    meta = next(
        (
            str(p.resolve())
            for p in (
                folder / "standardised.metadata.json",
                folder / f"{stem}.metadata.json",
                folder / "metadata.json",
            )
            if p.is_file()
        ),
        "",
    )
    qc = next(
        (
            str(p.resolve())
            for p in (
                folder / "standardised.qc.json",
                folder / f"{stem}.qc.json",
                folder / "qc.json",
            )
            if p.is_file()
        ),
        "",
    )
    return meta, qc


def preflight_row(
    raw: dict[str, Any],
    ld_population: str,
    allow_ancestry_mismatch: bool,
    min_maf: float = 0.01,
) -> tuple[dict[str, Any] | None, list[str]]:
    errors: list[str] = []
    vcf_str = (raw.get("vcf") or "").strip()
    if not vcf_str:
        return None, ["missing vcf path"]
    vcf = Path(vcf_str).expanduser().resolve()
    if not vcf.is_file():
        return None, [f"VCF file not found: {vcf}"]

    idx = find_index_for_vcf(vcf, raw.get("vcf_index") or None)
    if idx is None:
        return None, [f"Tabix/CSI index (.tbi/.csi) missing for {vcf}"]

    auto_meta, auto_qc = detect_sidecars_near_vcf(vcf)
    meta_str = (raw.get("metadata") or "").strip() or auto_meta
    qc_str = (raw.get("qc") or "").strip() or auto_qc

    meta: dict[str, Any] = {}
    qc: dict[str, Any] = {}
    if meta_str:
        mp = Path(meta_str).expanduser().resolve()
        if not mp.is_file():
            errors.append(f"metadata file not found: {mp}")
        else:
            try:
                meta = json.loads(mp.read_text(encoding="utf-8"))
                meta_str = str(mp)
            except Exception as exc:  # noqa: BLE001
                errors.append(f"invalid metadata JSON {mp}: {exc}")
    if qc_str:
        qp = Path(qc_str).expanduser().resolve()
        if not qp.is_file():
            errors.append(f"qc file not found: {qp}")
        else:
            try:
                qc = json.loads(qp.read_text(encoding="utf-8"))
                qc_str = str(qp)
            except Exception as exc:  # noqa: BLE001
                errors.append(f"invalid qc JSON {qp}: {exc}")

    if errors:
        return None, errors

    study_id = (
        (raw.get("study_id") or "").strip()
        or str(meta.get("study_id") or "").strip()
        or (vcf.parent.parent.parent.name if vcf.parent.parent.name == "standardise" else "study")
    )
    default_gwas = (
        vcf.parent.name
        if vcf.name.startswith("standardised.")
        else re.sub(r"\.(gwas\.)?vcf(\.gz|\.bgz)?$", "", vcf.name)
    )
    gwas_id = (
        (raw.get("gwas_id") or raw.get("output_id") or "").strip()
        or str(meta.get("output_id") or meta.get("gwas_id") or "").strip()
        or default_gwas
    )
    trait = (raw.get("trait") or "").strip() or str(meta.get("phenotype") or gwas_id)
    ancestry = (raw.get("ancestry") or "").strip() or str(meta.get("ancestry") or "")
    gene_id = (raw.get("gene_id") or raw.get("uniprot_id") or "").strip() or str(meta.get("uniprot_id") or "")
    row_pop = (raw.get("ld_pop") or raw.get("ld_population") or ld_population).strip()
    if row_pop not in REFERENCE_POPULATIONS:
        errors.append(f"unsupported LD population '{row_pop}'")

    for label, bval in (
        ("metadata.coordinate_build", meta.get("coordinate_build")),
        ("qc.coordinate_build", qc.get("coordinate_build")),
    ):
        if bval is not None and normalized_build(bval) != "GRCh38":
            errors.append(f"{label} is '{bval}' (expected GRCh38)")

    if (
        ancestry in REFERENCE_POPULATIONS
        and ancestry != row_pop
        and not allow_ancestry_mismatch
    ):
        errors.append(
            f"single-ancestry GWAS ({ancestry}) does not match LD population ({row_pop})"
        )

    # Lightweight VCF header + index check
    try:
        with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as gz:
            saw_chrom = False
            saw_fmt_af = False
            for line in gz:
                if line.startswith("##reference="):
                    ref_b = normalized_build(line.split("=", 1)[1].strip())
                    if ref_b and ref_b != "GRCh38":
                        errors.append(f"VCF ##reference is {ref_b} (expected GRCh38)")
                elif line.startswith("##FORMAT=<ID=AF,") or line.startswith("##FORMAT=<ID=AF>"):
                    saw_fmt_af = True
                elif line.startswith("#CHROM\t"):
                    saw_chrom = True
                    break
                elif not line.startswith("#"):
                    break
            if not saw_chrom:
                errors.append("VCF missing #CHROM header line")
            if float(min_maf or 0.0) > 0.0 and not saw_fmt_af:
                errors.append("VCF lacks FORMAT/AF required for --min-maf > 0 (pass --min-maf 0 to disable)")
    except OSError as exc:
        errors.append(f"cannot read gzipped VCF header: {exc}")

    res = subprocess.run(
        ["bcftools", "index", "--nrecords", str(vcf)],
        capture_output=True,
        text=True,
    )
    if res.returncode != 0:
        errors.append(f"bcftools index check failed: {res.stderr.strip()}")
    else:
        try:
            if int(res.stdout.strip()) < 1:
                errors.append("VCF index reports 0 records")
        except ValueError:
            errors.append("invalid record count from bcftools index")

    if errors:
        return {
            "study_id": study_id,
            "gwas_id": gwas_id,
            "vcf": str(vcf),
            "errors": errors,
        }, errors

    return {
        "study_id": study_id,
        "gwas_id": gwas_id,
        "row_key": f"{study_id}/{gwas_id}",
        "trait": trait,
        "ancestry": ancestry,
        "gene_id": gene_id,
        "ld_population": row_pop,
        "vcf": str(vcf),
        "vcf_index": str(idx),
        "metadata": meta_str,
        "qc": qc_str,
    }, []


def compute_settings_hash(row: dict[str, Any], params: dict[str, Any]) -> str:
    vcf_stat = Path(row["vcf"]).stat()
    payload = {
        "schema": 3,
        "vcf": row["vcf"],
        "vcf_size": vcf_stat.st_size,
        "ld_population": row["ld_population"],
        "gene_id": row["gene_id"],
        "methods": params["methods"],
        "output_format": params["output_format"],
        "significance_lp": round(float(params["significance_lp"]), 8),
        "cis_significance_lp": (
            round(float(params["cis_significance_lp"]), 8)
            if params.get("cis_significance_lp") is not None
            else None
        ),
        "region_flank_bp": int(params["region_flank_bp"]),
        "max_region_width_bp": int(params["max_region_width_bp"]),
        "gene_window_bp": int(params["gene_window_bp"]),
        "exclude_mhc": bool(params["exclude_mhc"]),
        "max_loci": int(params["max_loci"]),
        "max_causal": int(params["max_causal"]),
        "max_iter": int(params["max_iter"]),
        "cis_max_iter": int(params.get("cis_max_iter", 200)),
        "min_maf": round(float(params.get("min_maf", 0.01)), 8),
        "susie_refine": bool(params.get("susie_refine", False)),
        "effective_n": params.get("effective_n"),
        "reference_sample_size": params.get("reference_sample_size"),
    }
    raw = json.dumps(payload, sort_keys=True).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:16]


def load_state(state_file: Path) -> dict[str, dict[str, str]]:
    if not state_file.is_file():
        return {}
    with state_file.open(encoding="utf-8", newline="") as fh:
        return {row["row_key"]: row for row in csv.DictReader(fh, delimiter="\t") if row.get("row_key")}


def outputs_are_complete(outdir: Path, study_id: str, gwas_id: str,
                         settings_hash: str, state: dict[str, dict[str, str]],
                         require_report: bool, methods: list[str],
                         output_format: str, keep_intermediates: bool = False) -> bool:
    gdir = outdir / study_id / gwas_id
    if not (gdir / "qc.json").is_file() or not (gdir / "provenance.json").is_file():
        return False
    exts = ["rda", "rds"] if output_format == "both" else [output_format]
    for m in methods:
        stem = METHOD_STEMS[m]
        for ext in exts:
            if not (gdir / "methods" / f"{stem}.{ext}").is_file():
                return False
    if require_report and not (gdir / "report.html").is_file():
        return False
    rec = state.get(f"{study_id}/{gwas_id}")
    if rec and rec.get("settings_hash") == settings_hash and rec.get("status") in ("PASS", "WARN_DIAGNOSTIC_FLAGS", "NO_SIGNIFICANT_LOCI", "NO_PANEL_MATCHED_LOCI"):
        if keep_intermediates and int(rec.get("analyzed_regions", 0) or 0) > 0 and not (gdir / "loci").is_dir():
            return False
        return True
    return False


def atomic_write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as tmp:
        json.dump(data, tmp, indent=2, sort_keys=True)
        tmp.write("\n")
        tmp_path = Path(tmp.name)
    os.replace(tmp_path, path)


def cmd_prepare(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir).expanduser().resolve()
    state_file = Path(args.state_file).expanduser().resolve()
    rows_dir = Path(args.rows_dir).resolve()
    rows_dir.mkdir(parents=True, exist_ok=True)

    raw_rows: list[dict[str, Any]] = []
    source_desc = ""
    if args.vcf:
        source_desc = str(Path(args.vcf).expanduser().resolve())
        raw_rows.append({
            "vcf": args.vcf,
            "vcf_index": args.vcf_index or "",
            "metadata": args.metadata or "",
            "qc": args.qc or "",
            "study_id": args.study_id or "",
            "gwas_id": args.gwas_id or "",
            "trait": args.trait or "",
            "ancestry": args.ancestry or "",
            "gene_id": args.gene_id or "",
        })
    elif args.manifest:
        mf = Path(args.manifest).expanduser().resolve()
        source_desc = str(mf)
        with mf.open(encoding="utf-8", newline="") as fh:
            lines = [line for line in fh if line.strip() and not line.lstrip().startswith("#")]
        for r in csv.DictReader(lines, delimiter="\t"):
            if r.get("gwas_dir") and not r.get("vcf"):
                resolved = resolve_gwas_folder(Path(r["gwas_dir"]).expanduser().resolve(), r.get("study_id"))
                if resolved:
                    merged = {**resolved, **{k: v for k, v in r.items() if v}}
                    raw_rows.append(merged)
                else:
                    raw_rows.append(r)
            else:
                raw_rows.append(r)
    elif args.input_dir:
        source_desc = str(Path(args.input_dir).expanduser().resolve())
        raw_rows = discover_input_dir(Path(args.input_dir))
    else:
        raise ValueError("Supply one of --vcf, --manifest, or --input-dir.")

    rows_filtered = 0
    if args.only_study:
        filtered = []
        for r in raw_rows:
            sid = (r.get("study_id") or "").strip()
            if sid and sid != args.only_study:
                rows_filtered += 1
                continue
            filtered.append(r)
        raw_rows = filtered

    methods_list = normalize_methods_spec(args.methods)
    out_fmt = args.output_format.strip().lower()
    if out_fmt not in ("rda", "rds", "both"):
        raise ValueError(f"Unsupported --output-format '{out_fmt}'")

    if args.min_maf < 0.0 or args.min_maf >= 0.5:
        raise ValueError(f"--min-maf must be in [0.0, 0.5), got {args.min_maf}")

    params = {
        "methods": methods_list,
        "output_format": out_fmt,
        "significance_lp": args.significance_lp,
        "cis_significance_lp": args.cis_significance_lp,
        "min_maf": args.min_maf,
        "region_flank_bp": args.region_flank_bp,
        "max_region_width_bp": args.max_region_width_bp,
        "gene_window_bp": args.gene_window_bp,
        "exclude_mhc": args.exclude_mhc,
        "max_loci": args.max_loci,
        "max_causal": args.max_causal,
        "max_iter": args.max_iter,
        "cis_max_iter": args.cis_max_iter,
        "susie_refine": args.susie_refine,
        "effective_n": args.effective_n,
        "reference_sample_size": args.reference_sample_size,
    }

    state = load_state(state_file)
    seen_keys: set[str] = set()
    eligible: list[dict[str, Any]] = []
    invalid: list[dict[str, Any]] = []
    skipped_complete = 0

    for raw in raw_rows:
        checked, errs = preflight_row(
            raw=raw,
            ld_population=args.population,
            allow_ancestry_mismatch=args.allow_ancestry_mismatch,
            min_maf=args.min_maf,
        )
        if errs:
            invalid.append({
                "study_id": (checked or raw).get("study_id", ""),
                "gwas_id": (checked or raw).get("gwas_id", ""),
                "vcf": raw.get("vcf", ""),
                "errors": errs,
            })
            continue
        assert checked is not None
        if args.only_study and checked["study_id"] != args.only_study:
            rows_filtered += 1
            continue
        if checked["row_key"] in seen_keys:
            invalid.append({
                "study_id": checked["study_id"],
                "gwas_id": checked["gwas_id"],
                "vcf": checked["vcf"],
                "errors": [f"duplicate study_id/gwas_id '{checked['row_key']}' in run"],
            })
            continue
        seen_keys.add(checked["row_key"])

        shash = compute_settings_hash(checked, params)
        checked["settings_hash"] = shash
        if not args.force and outputs_are_complete(
            outdir, checked["study_id"], checked["gwas_id"], shash, state,
            args.report, methods_list, out_fmt, args.keep_intermediates
        ):
            skipped_complete += 1
            continue
        eligible.append(checked)

    deferred = 0
    if args.row_limit > 0 and len(eligible) > args.row_limit:
        deferred = len(eligible) - args.row_limit
        to_run = eligible[: args.row_limit]
    else:
        to_run = eligible

    manifest_out = Path(args.output).resolve()
    manifest_out.parent.mkdir(parents=True, exist_ok=True)
    with manifest_out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RUN_MANIFEST_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for idx, row in enumerate(to_run, start=1):
            safe_slug = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{row['study_id']}__{row['gwas_id']}")[:80]
            row_id = f"row_{idx:05d}_{safe_slug}"
            row["row_id"] = row_id
            writer.writerow(row)
            if not args.dry_run:
                row_gwas_dir = outdir / row["study_id"] / row["gwas_id"]
                if args.force and (row_gwas_dir / "loci").is_dir():
                    shutil.rmtree(row_gwas_dir / "loci", ignore_errors=True)
                row_log_dir = row_gwas_dir / "logs"
                if row_log_dir.is_dir():
                    stale_prefixes = ("01_", "02_", "03_", "04_") if args.force else ("01_", "04_")
                    for stale in row_log_dir.iterdir():
                        if stale.is_file() and (
                            stale.name.startswith(stale_prefixes)
                            or stale.name == "trace.tsv"
                        ):
                            stale.unlink()
                payload = {
                    **row,
                    "outdir": str(outdir),
                    "chain": str(Path(args.chain).expanduser().resolve()),
                    "grch37_fasta": str(Path(args.grch37_fasta).expanduser().resolve()),
                    "gene_annotation": (
                        str(Path(args.gene_annotation).expanduser().resolve())
                        if args.gene_annotation and Path(args.gene_annotation).expanduser().is_file()
                        else ""
                    ),
                    "panukb_dir": (
                        str(Path(args.panukb_dir).expanduser().resolve())
                        if args.panukb_dir
                        else ""
                    ),
                    "polyfun_dir": (
                        str(Path(args.polyfun_dir).expanduser().resolve())
                        if args.polyfun_dir
                        else ""
                    ),
                    "local_ld_panel": (
                        str(Path(args.local_ld_panel).expanduser().resolve())
                        if args.local_ld_panel
                        else ""
                    ),
                    "bm_path": (
                        str(Path(args.bm_path).expanduser().resolve())
                        if args.bm_path
                        else ""
                    ),
                    "cache_dir": str(Path(args.cache_dir).expanduser().resolve()),
                    "force": args.force,
                    "report": args.report,
                    "keep_intermediates": args.keep_intermediates,
                    "allow_ancestry_mismatch": args.allow_ancestry_mismatch,
                    **params,
                }
                atomic_write_json(rows_dir / f"{row_id}.json", payload)

    summary = {
        "source": source_desc,
        "outdir": str(outdir),
        "state_file": str(state_file),
        "dry_run": args.dry_run,
        "rows_discovered": len(raw_rows) + rows_filtered,
        "rows_filtered": rows_filtered,
        "rows_skipped_complete": skipped_complete,
        "rows_eligible": len(eligible),
        "rows_deferred_due_to_limit": deferred,
        "rows_to_run": len(to_run),
        "invalid": invalid,
    }
    atomic_write_json(Path(args.summary).resolve(), summary)
    return 0


def cmd_update_state(args: argparse.Namespace) -> int:
    state_file = Path(args.state_file).expanduser().resolve()
    state_file.parent.mkdir(parents=True, exist_ok=True)
    lock_path = state_file.parent / f"{state_file.name}.lock"

    updated = 0
    failed = 0
    completed_keys: list[str] = []
    failed_rows: list[dict[str, Any]] = []

    with open(lock_path, "a+", encoding="utf-8") as lock_fh:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
        state = load_state(state_file)

        for rec_path_str in args.records:
            rec_path = Path(rec_path_str)
            if not rec_path.is_file():
                continue
            rec = json.loads(rec_path.read_text(encoding="utf-8"))
            rkey = rec["row_key"]
            if rec.get("ok"):
                state[rkey] = {
                    "row_key": rkey,
                    "study_id": rec["study_id"],
                    "gwas_id": rec["gwas_id"],
                    "settings_hash": rec["settings_hash"],
                    "status": rec.get("overall_status", "PASS"),
                    "analyzed_regions": str(rec.get("analyzed_regions", 0)),
                    "cs_rss": str(rec.get("cs_rss", 0)),
                    "cs_rrss": str(rec.get("cs_rrss", 0)),
                    "output_dir": rec.get("output_dir", ""),
                }
                updated += 1
                completed_keys.append(rkey)
            else:
                failed += 1
                failed_rows.append({
                    "row_key": rkey,
                    "study_id": rec.get("study_id", ""),
                    "gwas_id": rec.get("gwas_id", ""),
                    "error": rec.get("error", "task failed"),
                })

        if updated > 0 or state_file.is_file():
            with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=state_file.parent, newline="", delete=False) as tmp:
                writer = csv.DictWriter(tmp, fieldnames=STATE_COLUMNS, delimiter="\t")
                writer.writeheader()
                for k in sorted(state):
                    writer.writerow(state[k])
                tmp_path = Path(tmp.name)
            os.replace(tmp_path, state_file)

    if args.outdir and not getattr(args, "skip_manifest_report", False):
        manifest_script = Path(__file__).with_name("render_manifest_report.py")
        if manifest_script.is_file():
            subprocess.run(
                [sys.executable, str(manifest_script), "--outdir", str(Path(args.outdir).resolve())],
                check=False,
            )

    if args.output:
        atomic_write_json(
            Path(args.output).resolve(),
            {
                "updated": updated,
                "failed": failed,
                "completed_keys": completed_keys,
                "failed_rows": failed_rows,
            },
        )
    return 0


def cmd_copy_logs(args: argparse.Namespace) -> int:
    outdir = Path(args.outdir).expanduser().resolve()
    if args.row_json and Path(args.row_json).is_file():
        row = json.loads(Path(args.row_json).read_text(encoding="utf-8"))
        log_dir = outdir / row["study_id"] / row["gwas_id"] / "logs"
        step_stem = re.sub(r"[^A-Za-z0-9._-]+", "_", (args.step_name or row["gwas_id"]).strip())
        try:
            log_dir.mkdir(parents=True, exist_ok=True)
            wrote_log = False
            for cand_name in (".command.log", ".command.out"):
                p = Path(cand_name)
                if p.is_file() and p.stat().st_size > 0:
                    shutil.copy2(p, log_dir / f"{step_stem}.log")
                    wrote_log = True
                    break
            if not wrote_log:
                (log_dir / f"{step_stem}.log").write_text(f"[completed] {step_stem}\n", encoding="utf-8")
            err_p = Path(".command.err")
            if err_p.is_file() and err_p.stat().st_size > 0:
                shutil.copy2(err_p, log_dir / f"{step_stem}.err")
        except OSError as exc:
            sys.stderr.write(f"[WARN] Could not write step log to {log_dir}: {exc}\n")
        return 0

    manifest_rows: list[dict[str, str]] = []
    if args.run_manifest and Path(args.run_manifest).is_file():
        with Path(args.run_manifest).open(encoding="utf-8", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("study_id") and r.get("gwas_id"):
                    manifest_rows.append(r)

    trace_header = ""
    trace_lines: list[str] = []
    if args.trace_file and Path(args.trace_file).is_file():
        raw_trace = Path(args.trace_file).read_text(encoding="utf-8", errors="replace").splitlines()
        if raw_trace:
            trace_header = raw_trace[0]
            trace_lines = raw_trace[1:]

    for r in manifest_rows:
        log_dir = outdir / r["study_id"] / r["gwas_id"] / "logs"
        try:
            log_dir.mkdir(parents=True, exist_ok=True)
            if args.nextflow_log and Path(args.nextflow_log).is_file():
                shutil.copy2(args.nextflow_log, log_dir / Path(args.nextflow_log).name)
            if args.slurm_out and Path(args.slurm_out).is_file():
                shutil.copy2(args.slurm_out, log_dir / Path(args.slurm_out).name)
            row_id = r.get("row_id", "")
            if trace_header and row_id:
                matched_trace = [ln for ln in trace_lines if row_id in ln]
                if matched_trace:
                    (log_dir / "trace.tsv").write_text(
                        "\n".join([trace_header, *matched_trace]) + "\n",
                        encoding="utf-8",
                    )
        except OSError as exc:
            sys.stderr.write(f"[WARN] Could not write run log to {log_dir}: {exc}\n")
    return 0


def str2bool(v: str) -> bool:
    return str(v).strip().lower() in ("1", "true", "yes", "y")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_prep = sub.add_parser("prepare")
    p_prep.add_argument("--input-dir")
    p_prep.add_argument("--manifest")
    p_prep.add_argument("--vcf")
    p_prep.add_argument("--vcf-index")
    p_prep.add_argument("--metadata")
    p_prep.add_argument("--qc")
    p_prep.add_argument("--study-id")
    p_prep.add_argument("--gwas-id")
    p_prep.add_argument("--trait")
    p_prep.add_argument("--ancestry")
    p_prep.add_argument("--gene-id")
    p_prep.add_argument("--only-study")
    p_prep.add_argument("--outdir", required=True)
    p_prep.add_argument("--state-file", required=True)
    p_prep.add_argument("--output", required=True)
    p_prep.add_argument("--summary", required=True)
    p_prep.add_argument("--rows-dir", required=True)
    p_prep.add_argument("--chain", required=True)
    p_prep.add_argument("--grch37-fasta", required=True)
    p_prep.add_argument("--gene-annotation")
    p_prep.add_argument("--panukb-dir")
    p_prep.add_argument("--polyfun-dir")
    p_prep.add_argument("--local-ld-panel")
    p_prep.add_argument("--bm-path")
    p_prep.add_argument("--cache-dir", required=True)
    p_prep.add_argument("--population", default="EUR")
    p_prep.add_argument("--methods", default="p_ld,susie_rss,susie_rrss,carma,polyfun")
    p_prep.add_argument("--output-format", default="rda")
    p_prep.add_argument("--reference-sample-size", type=int)
    p_prep.add_argument("--effective-n", type=float)
    p_prep.add_argument("--significance-lp", type=float, default=7.301029995664)
    p_prep.add_argument("--cis-significance-lp", type=float)
    p_prep.add_argument("--min-maf", type=float, default=0.01)
    p_prep.add_argument("--region-flank-bp", type=int, default=1_500_000)
    p_prep.add_argument("--max-region-width-bp", type=int, default=9_500_000)
    p_prep.add_argument("--gene-window-bp", type=int, default=1_000_000)
    p_prep.add_argument("--exclude-mhc", type=str2bool, default=False)
    p_prep.add_argument("--max-loci", type=int, default=0)
    p_prep.add_argument("--max-causal", type=int, default=10)
    p_prep.add_argument("--max-iter", type=int, default=100)
    p_prep.add_argument("--cis-max-iter", type=int, default=200)
    p_prep.add_argument("--susie-refine", type=str2bool, default=False)
    p_prep.add_argument("--row-limit", type=int, default=0)
    p_prep.add_argument("--force", type=str2bool, default=False)
    p_prep.add_argument("--dry-run", type=str2bool, default=False)
    p_prep.add_argument("--report", type=str2bool, default=True)
    p_prep.add_argument("--keep-intermediates", type=str2bool, default=False)
    p_prep.add_argument("--allow-ancestry-mismatch", type=str2bool, default=False)

    p_upd = sub.add_parser("update-state")
    p_upd.add_argument("--state-file", required=True)
    p_upd.add_argument("--outdir")
    p_upd.add_argument("--output", default="")
    p_upd.add_argument("--skip-manifest-report", action="store_true")
    p_upd.add_argument("--records", nargs="*", default=[])

    p_log = sub.add_parser("copy-logs")
    p_log.add_argument("--outdir", required=True)
    p_log.add_argument("--row-json")
    p_log.add_argument("--step-name")
    p_log.add_argument("--run-manifest")
    p_log.add_argument("--nextflow-log")
    p_log.add_argument("--trace-file")
    p_log.add_argument("--slurm-out")

    args = parser.parse_args()
    if args.cmd == "prepare":
        return cmd_prepare(args)
    if args.cmd == "update-state":
        return cmd_update_state(args)
    return cmd_copy_logs(args)


if __name__ == "__main__":
    raise SystemExit(main())
