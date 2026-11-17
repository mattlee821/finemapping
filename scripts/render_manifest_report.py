#!/usr/bin/env python3
"""Render a self-contained batch manifest.html index across fine-mapped GWAS outputs."""

from __future__ import annotations

import argparse
import html
import json
import os
from pathlib import Path
from typing import Any


def badge(label: str, tone: str = "neutral") -> str:
    return f'<span class="badge badge-{tone}">{html.escape(label)}</span>'


def collect_runs(outdir: Path) -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    for qc_file in sorted(outdir.rglob("qc.json")):
        gwas_dir = qc_file.parent
        prov_file = gwas_dir / "provenance.json"
        meta_file = gwas_dir / "metadata.json"
        if not (prov_file.is_file() and meta_file.is_file()):
            continue
        qc = json.loads(qc_file.read_text(encoding="utf-8"))
        prov = json.loads(prov_file.read_text(encoding="utf-8"))
        meta = json.loads(meta_file.read_text(encoding="utf-8"))
        if "analyzed_regions_count" not in qc:
            continue
        report_file = gwas_dir / "report.html"
        rel_report = (
            os.path.relpath(report_file, outdir)
            if report_file.is_file()
            else None
        )
        rda_bytes = sum(
            f.stat().st_size
            for f in (gwas_dir / "methods").glob("*.rda")
            if f.is_file()
        )
        runs.append({
            "study_id": prov.get("study_id", ""),
            "gwas_id": prov.get("gwas_id", ""),
            "trait": prov.get("trait") or meta.get("phenotype") or prov.get("gwas_id", ""),
            "ld_pop": prov.get("ld_population", ""),
            "status": qc.get("overall_status", "UNKNOWN"),
            "sig_variants": qc.get("genomewide_significant_variants", 0),
            "regions": qc.get("analyzed_regions_count", 0),
            "cs_rss": qc.get("total_credible_sets_susie_rss", 0),
            "cs_rrss": qc.get("total_credible_sets_susie_rrss", 0),
            "cis_status": (qc.get("cis_gene") or {}).get("status", "not_applicable"),
            "mhc_status": (qc.get("mhc") or {}).get("status", "no_significant_variants"),
            "nonconverged": len(qc.get("nonconverged_regions") or []),
            "mismatch": len(qc.get("mismatch_flagged_regions") or []),
            "runtime_s": prov.get("runtime_seconds", 0.0),
            "rda_mb": rda_bytes / (1024 * 1024),
            "rel_report": rel_report,
        })
    return runs


def render_manifest_html(outdir: Path) -> Path:
    runs = collect_runs(outdir)
    n_pass = sum(1 for r in runs if r["status"] == "PASS")
    n_warn = sum(1 for r in runs if r["status"] == "WARN_DIAGNOSTIC_FLAGS")
    n_null = sum(1 for r in runs if r["status"] == "NO_SIGNIFICANT_LOCI")
    total_mb = sum(r["rda_mb"] for r in runs)

    rows_html = []
    for r in runs:
        tone = {
            "PASS": "ok",
            "WARN_DIAGNOSTIC_FLAGS": "warn",
            "NO_SIGNIFICANT_LOCI": "neutral",
        }.get(r["status"], "bad")
        link = (
            f'<a href="{html.escape(r["rel_report"])}"><code>{html.escape(r["gwas_id"])}</code></a>'
            if r["rel_report"]
            else f'<code>{html.escape(r["gwas_id"])}</code>'
        )
        rows_html.append(
            "<tr>"
            f"<td><code>{html.escape(str(r['study_id']))}</code></td>"
            f"<td>{link}</td>"
            f"<td>{html.escape(str(r['trait']))}</td>"
            f"<td><code>{html.escape(str(r['ld_pop']))}</code></td>"
            f"<td>{badge(str(r['status']), tone)}</td>"
            f"<td class='num'>{r['sig_variants']:,}</td>"
            f"<td class='num'>{r['regions']:,}</td>"
            f"<td class='num'>{r['cs_rss']:,} / {r['cs_rrss']:,}</td>"
            f"<td><code>{html.escape(str(r['cis_status']))}</code></td>"
            f"<td><code>{html.escape(str(r['mhc_status']))}</code></td>"
            f"<td class='num'>{r['nonconverged']} / {r['mismatch']}</td>"
            f"<td class='num'>{r['rda_mb']:.2f} MB</td>"
            f"<td class='num'>{r['runtime_s']:.1f}s</td>"
            "</tr>"
        )

    body_rows = "\n".join(rows_html) if rows_html else (
        '<tr><td colspan="13" class="muted">No fine-mapped GWAS outputs found.</td></tr>'
    )

    doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>GWASfinemapping Manifest Summary</title>
<style>
:root {{
  --bg: #f8fafc; --card: #ffffff; --fg: #0f172a; --muted: #64748b;
  --border: #e2e8f0; --link: #0f766e;
  --ok-bg: #dcfce7; --ok-fg: #166534;
  --warn-bg: #fef3c7; --warn-fg: #92400e;
  --bad-bg: #fee2e2; --bad-fg: #991b1b;
  --neu-bg: #e2e8f0; --neu-fg: #334155;
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --bg: #0f172a; --card: #1e293b; --fg: #f1f5f9; --muted: #94a3b8;
    --border: #334155; --link: #2dd4bf;
    --ok-bg: #14532d; --ok-fg: #bbf7d0;
    --warn-bg: #78350f; --warn-fg: #fde68a;
    --bad-bg: #7f1d1d; --bad-fg: #fecaca;
    --neu-bg: #334155; --neu-fg: #e2e8f0;
  }}
}}
body {{
  margin: 0; padding: 24px; background: var(--bg); color: var(--fg);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  font-size: 14px;
}}
.container {{ max-width: 1400px; margin: 0 auto; }}
a {{ color: var(--link); text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px; margin: 16px 0 22px; }}
.card, .panel {{ background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 14px 16px; }}
.label {{ font-size: 11.5px; text-transform: uppercase; color: var(--muted); }}
.val {{ font-size: 20px; font-weight: 600; margin-top: 4px; }}
.muted {{ color: var(--muted); }}
.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.badge {{ display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11.5px; font-weight: 600; }}
.badge-ok {{ background: var(--ok-bg); color: var(--ok-fg); }}
.badge-warn {{ background: var(--warn-bg); color: var(--warn-fg); }}
.badge-bad {{ background: var(--bad-bg); color: var(--bad-fg); }}
.badge-neutral {{ background: var(--neu-bg); color: var(--neu-fg); }}
.table-wrap {{ overflow-x: auto; }}
table {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
th, td {{ padding: 8px 10px; border-bottom: 1px solid var(--border); text-align: left; }}
th {{ color: var(--muted); font-size: 12px; background: rgba(148, 163, 184, 0.08); }}
</style>
</head>
<body>
<div class="container">
  <h1>GWASfinemapping Manifest Summary</h1>
  <div class="muted">Output directory: <code>{html.escape(str(outdir))}</code></div>
  <div class="grid">
    <div class="card"><div class="label">Total GWAS</div><div class="val">{len(runs):,}</div></div>
    <div class="card"><div class="label">PASS / WARN / No Loci</div><div class="val">{n_pass} / {n_warn} / {n_null}</div></div>
    <div class="card"><div class="label">Total .rda Footprint</div><div class="val">{total_mb:.2f} MB</div></div>
  </div>
  <div class="panel table-wrap">
    <table>
      <thead>
        <tr>
          <th>Study</th><th>GWAS (Report)</th><th>Trait</th><th>LD Pop</th><th>Status</th>
          <th class="num">Sig SNPs</th><th class="num">Regions</th><th class="num">CS (RSS / rRSS)</th>
          <th>cis Gate</th><th>MHC</th><th class="num">Nonconv / Mismatch</th>
          <th class="num">.rda Size</th><th class="num">Time</th>
        </tr>
      </thead>
      <tbody>{body_rows}</tbody>
    </table>
  </div>
</div>
</body>
</html>
"""
    target = outdir / "manifest.html"
    target.write_text(doc, encoding="utf-8")
    return target


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, help="Root output directory containing <study>/<gwas>/ folders")
    args = parser.parse_args()
    out = render_manifest_html(Path(args.outdir).resolve())
    print(out)


if __name__ == "__main__":
    main()
