nextflow.enable.dsl = 2

def opt(String flag, def val) {
    val != null && val.toString() != '' ? "${flag} '${val}'" : ''
}

def parseMethodsList(String spec) {
    def raw = (spec ?: 'p_ld,susie_rss,susie_rrss,carma,polyfun').split(/[\s,]+/).collect { token ->
        def m = token.trim().toLowerCase().replace('-', '_')
        m == 'polyfun_susie' ? 'polyfun' : m
    }.findAll { it }
    return raw.unique()
}

process PREPARE_MANIFEST {
    tag "prepare-manifest"
    cache false
    publishDir "${params.log_root ?: params.outdir}", mode: 'copy', overwrite: true, pattern: 'run_manifest*'

    output:
    path 'run_manifest.tsv', emit: manifest
    path 'run_manifest_summary.json', emit: summary
    path 'rows/*.json', emit: rows, optional: true

    script:
    """
    python3 '${projectDir}/scripts/prepare_manifest.py' prepare \\
      ${opt('--input-dir', params.input_dir)} \\
      ${opt('--manifest', params.manifest)} \\
      ${opt('--vcf', params.vcf)} \\
      ${opt('--vcf-index', params.vcf_index)} \\
      ${opt('--metadata', params.metadata)} \\
      ${opt('--qc', params.qc)} \\
      ${opt('--study-id', params.study_id)} \\
      ${opt('--gwas-id', params.gwas_id)} \\
      ${opt('--trait', params.trait)} \\
      ${opt('--ancestry', params.ancestry)} \\
      ${opt('--gene-id', params.gene_id)} \\
      ${opt('--only-study', params.only_study)} \\
      --outdir '${params.outdir}' \\
      --state-file '${params.outdir}/.finemap_state.tsv' \\
      --output run_manifest.tsv \\
      --summary run_manifest_summary.json \\
      --rows-dir rows \\
      --chain '${params.chain}' \\
      --grch37-fasta '${params.grch37_fasta}' \\
      ${opt('--gene-annotation', params.gene_annotation)} \\
      ${opt('--panukb-dir', params.panukb_dir)} \\
      ${opt('--polyfun-dir', params.polyfun_dir)} \\
      ${opt('--local-ld-panel', params.local_ld_panel)} \\
      ${opt('--bm-path', params.bm_path)} \\
      --cache-dir '${params.cache_dir}' \\
      --population '${params.population}' \\
      --methods '${params.methods}' \\
      --output-format '${params.output_format}' \\
      ${opt('--reference-sample-size', params.reference_sample_size)} \\
      ${opt('--effective-n', params.effective_n)} \\
      --significance-lp '${params.significance_lp}' \\
      ${opt('--cis-significance-lp', params.cis_significance_lp)} \\
      --region-flank-bp '${params.region_flank_bp}' \\
      --max-region-width-bp '${params.max_region_width_bp}' \\
      --gene-window-bp '${params.gene_window_bp}' \\
      --exclude-mhc '${params.exclude_mhc}' \\
      --max-loci '${params.max_loci}' \\
      --max-causal '${params.max_causal}' \\
      --max-iter '${params.max_iter}' \\
      --row-limit '${params.row_limit}' \\
      --force '${params.force}' \\
      --dry-run '${params.dry_run}' \\
      --report '${params.report}' \\
      --keep-intermediates '${params.keep_intermediates}' \\
      --allow-ancestry-mismatch '${params.allow_ancestry_mismatch}'
    """
}

process SCREEN_GWAS {
    tag "${row_id}"
    cache false

    input:
    tuple val(row_id), path(row_json)

    output:
    tuple val(row_id), path(row_json), path('screen_summary.json'), path('loci_list.tsv')

    script:
    """
    python3 - '${projectDir}/scripts/screen_gwas.py' '${row_json}' <<'PY'
import csv, json, subprocess, sys
from pathlib import Path

screen_script = sys.argv[1]
row = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))

cmd = [
    sys.executable, screen_script,
    "--vcf", row["vcf"],
    "--study-id", row["study_id"],
    "--gwas-id", row["gwas_id"],
    "--chain", row["chain"],
    "--population", row["ld_population"],
    "--significance-lp", str(row["significance_lp"]),
    "--region-flank-bp", str(row["region_flank_bp"]),
    "--gene-window-bp", str(row["gene_window_bp"]),
    "--max-region-width-bp", str(row["max_region_width_bp"]),
    "--max-loci", str(row["max_loci"]),
    "--out", "screen_summary.json",
]
if row.get("metadata"):
    cmd.extend(["--metadata", row["metadata"]])
if row.get("qc"):
    cmd.extend(["--qc", row["qc"]])
if row.get("trait"):
    cmd.extend(["--trait", row["trait"]])
if row.get("ancestry"):
    cmd.extend(["--ancestry", row["ancestry"]])
if row.get("gene_id"):
    cmd.extend(["--gene-id", row["gene_id"]])
if row.get("gene_annotation"):
    cmd.extend(["--gene-annotation", row["gene_annotation"]])
if row.get("reference_sample_size") is not None:
    cmd.extend(["--reference-sample-size", str(row["reference_sample_size"])])
if row.get("cis_significance_lp") is not None:
    cmd.extend(["--cis-significance-lp", str(row["cis_significance_lp"])])
if row.get("exclude_mhc"):
    cmd.append("--exclude-mhc")
if row.get("allow_ancestry_mismatch"):
    cmd.append("--allow-ancestry-mismatch")

subprocess.run(cmd, check=True)

summary = json.loads(Path("screen_summary.json").read_text(encoding="utf-8"))
with open("loci_list.tsv", "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\\t")
    w.writerow(["locus_index", "region_id", "region_type"])
    for idx, reg in enumerate(summary.get("selected_regions") or []):
        w.writerow([idx, reg["region_id"], reg["region_type"]])
PY
    """
}

process PREP_LOCUS {
    tag "${row_id}:${region_id}"
    label 'stage_prep'
    cache false

    input:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(locus_index), val(region_id), val(region_type)

    output:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(region_id), val(region_type), env('SCRATCH_LOCI_ROOT')

    script:
    def scratch_base = "${params.cache_dir}/scratch/${workflow.sessionId}/${row_id}/loci"
    """
    python3 - '${projectDir}/scripts/prep_locus.py' '${row_json}' '${screen_json}' '${locus_index}' '${region_id}' '${scratch_base}' <<'PY'
import json, subprocess, sys
from pathlib import Path

prep_script, row_p, screen_p, loc_idx, reg_id, scratch_base = sys.argv[1:7]
row = json.loads(Path(row_p).read_text(encoding="utf-8"))
if row.get("keep_intermediates"):
    loci_root = Path(row["outdir"]) / row["study_id"] / row["gwas_id"] / "loci"
else:
    loci_root = Path(scratch_base)

outdir = loci_root / reg_id
cmd = [
    sys.executable, prep_script,
    "--screen-json", screen_p,
    "--locus-index", str(loc_idx),
    "--outdir", str(outdir),
    "--chain", row["chain"],
    "--grch37-fasta", row["grch37_fasta"],
    "--cache-dir", row["cache_dir"],
]
if row.get("panukb_dir"):
    cmd.extend(["--panukb-dir", row["panukb_dir"]])
if row.get("local_ld_panel"):
    cmd.extend(["--local-ld-panel", row["local_ld_panel"]])
if row.get("bm_path"):
    cmd.extend(["--bm-path", row["bm_path"]])
if row.get("keep_intermediates"):
    cmd.append("--keep-intermediates")

subprocess.run(cmd, check=True)
Path(".scratch_loci_root").write_text(str(loci_root), encoding="utf-8")
PY
    SCRATCH_LOCI_ROOT="\$(cat .scratch_loci_root)"
    """
}

process RUN_P_LD {
    tag "${row_id}:${region_id}:p_ld"
    label 'method_light'
    cache false

    input:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(region_id), val(region_type), val(scratch_loci_root)

    output:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(region_id), val(scratch_loci_root)

    script:
    def loc_dir = "${scratch_loci_root}/${region_id}"
    """
    Rscript '${projectDir}/scripts/methods/run_p_ld.R' \\
      --variants '${loc_dir}/matched_variants.tsv' \\
      --matrix '${loc_dir}/matrix.npy' \\
      --region-id '${region_id}' \\
      --region-type '${region_type}' \\
      --out-rds '${loc_dir}/methods/p_ld/region_result.rds' \\
      --qc-json '${loc_dir}/methods/p_ld/method_qc.json' \\
      --allow-skipped-locus true
    """
}

process RUN_BAYESIAN_METHOD {
    tag "${row_id}:${region_id}:${method_id}"
    label 'method_heavy'
    cache false

    input:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(region_id), val(region_type), val(scratch_loci_root), val(method_id)

    output:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(region_id), val(scratch_loci_root)

    script:
    def loc_dir = "${scratch_loci_root}/${region_id}"
    """
    python3 - '${projectDir}' '${row_json}' '${screen_json}' '${loc_dir}' '${region_id}' '${region_type}' '${method_id}' <<'PY'
import json, subprocess, sys
from pathlib import Path

proj = Path(sys.argv[1])
row = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
screen = json.loads(Path(sys.argv[3]).read_text(encoding="utf-8"))
loc_dir = Path(sys.argv[4])
region_id, region_type, method_id = sys.argv[5:8]

script_path = proj / "scripts" / "methods" / f"run_{method_id}.R"
if not script_path.is_file():
    raise FileNotFoundError(f"Standalone method script not found: {script_path}")

cmd = [
    "Rscript", str(script_path),
    "--variants", str(loc_dir / "matched_variants.tsv"),
    "--matrix", str(loc_dir / "matrix.npy"),
    "--region-id", region_id,
    "--region-type", region_type,
    "--out-rds", str(loc_dir / "methods" / method_id / "region_result.rds"),
    "--qc-json", str(loc_dir / "methods" / method_id / "method_qc.json"),
    "--max-causal", str(row.get("max_causal", 10)),
    "--max-iter", str(row.get("max_iter", 50)),
    "--reference-sample-size", str(screen["reference_sample_size"]),
    "--allow-skipped-locus", "true",
]
if row.get("effective_n") is not None:
    cmd.extend(["--effective-n", str(row["effective_n"])])
if method_id == "polyfun" and row.get("polyfun_dir"):
    cmd.extend(["--polyfun-priors", str(row["polyfun_dir"])])

subprocess.run(cmd, check=True)
PY
    """
}

process CLEAN_LOCUS_SCRATCH {
    tag "${row_id}:${region_id}"
    label 'method_light'
    cache false

    input:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(region_id), val(scratch_loci_root)

    output:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(scratch_loci_root)

    script:
    def loc_dir = "${scratch_loci_root}/${region_id}"
    """
    python3 - '${row_json}' '${loc_dir}' <<'PY'
import json, sys
from pathlib import Path

row = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
loc_dir = Path(sys.argv[2])
if not row.get("keep_intermediates"):
    for fn in ("matrix.npy", "matrix.ids.txt"):
        p = loc_dir / fn
        if p.exists():
            p.unlink()
PY
    """
}

process FINALISE_GWAS {
    tag "${row_id}"
    cache false
    afterScript "python3 '${projectDir}/scripts/prepare_manifest.py' copy-logs --outdir '${params.outdir}' --row-json '${row_json}' || true"

    input:
    tuple val(row_id), path(row_json), path(screen_json), val(scratch_loci_root)

    output:
    path '*.row_record.json', emit: records

    script:
    """
    python3 - '${projectDir}/scripts/finalise_gwas.py' '${row_id}' '${row_json}' '${screen_json}' '${scratch_loci_root}' '${params.compress}' <<'PY'
import json, subprocess, sys
from pathlib import Path

fin_script, row_id, row_p, screen_p, loci_root, compress = sys.argv[1:7]
row = json.loads(Path(row_p).read_text(encoding="utf-8"))
methods_csv = ",".join(row["methods"])

cmd = [
    sys.executable, fin_script,
    "--screen-json", screen_p,
    "--output-root", row["outdir"],
    "--methods", methods_csv,
    "--output-format", row.get("output_format", "rda"),
    "--compress", compress,
    "--chain", row["chain"],
    "--grch37-fasta", row["grch37_fasta"],
    "--max-causal", str(row.get("max_causal", 10)),
    "--max-iter", str(row.get("max_iter", 50)),
    "--report", "true" if row.get("report", True) else "false",
    "--settings-hash", row.get("settings_hash", ""),
    "--output-record", f"{row_id}.row_record.json",
]
if loci_root and loci_root != "NONE":
    cmd.extend(["--loci-dir", loci_root])
if row.get("effective_n") is not None:
    cmd.extend(["--effective-n", str(row["effective_n"])])
if row.get("keep_intermediates"):
    cmd.append("--keep-intermediates")
else:
    cmd.append("--cleanup-loci-scratch")
if row.get("force"):
    cmd.append("--force")

try:
    subprocess.run(cmd, check=True)
except Exception as exc:
    Path(f"{row_id}.row_record.json").write_text(json.dumps({
        "ok": False,
        "row_key": row["row_key"],
        "study_id": row["study_id"],
        "gwas_id": row["gwas_id"],
        "settings_hash": row.get("settings_hash", ""),
        "error": str(exc),
    }, indent=2) + "\\n", encoding="utf-8")
PY
    """
}

process UPDATE_STATE {
    tag "update-state"
    cache false
    publishDir "${params.log_root ?: params.outdir}", mode: 'copy', overwrite: true

    input:
    path records

    output:
    path 'state_update_summary.json', emit: summary

    script:
    def rec_args = records ? "--records ${records}" : ""
    """
    python3 '${projectDir}/scripts/prepare_manifest.py' update-state \\
      --state-file '${params.outdir}/.finemap_state.tsv' \\
      --outdir '${params.outdir}' \\
      --output state_update_summary.json \\
      ${rec_args}
    """
}

workflow {
    def ref_required = ['chain', 'grch37_fasta', 'grch37_fai', 'outdir']
    def missing_ref = ref_required.findAll { !params[it] }
    if (missing_ref) {
        error "Missing required parameter(s): ${missing_ref.collect { "--${it}" }.join(', ')}"
    }
    if (!params.input_dir && !params.manifest && !params.vcf) {
        error "Supply one of --input_dir <dir>, --manifest <tsv>, or --vcf <path>"
    }

    def methods_list = parseMethodsList(params.methods as String)

    PREPARE_MANIFEST()
    row_files_ch = PREPARE_MANIFEST.out.rows.flatten().map { rj -> tuple(rj.baseName, rj) }

    SCREEN_GWAS(row_files_ch)

    branched_ch = SCREEN_GWAS.out.branch { row_id, row_json, screen_json, loci_tsv ->
        def lines = loci_tsv.readLines().drop(1).findAll { it.trim() }
        has_loci: lines.size() > 0
            return tuple(row_id, row_json, screen_json, lines)
        zero_loci: true
            return tuple(row_id, row_json, screen_json, 'NONE')
    }

    locus_tasks_ch = branched_ch.has_loci.flatMap { row_id, row_json, screen_json, lines ->
        def n_loci = lines.size()
        lines.collect { line ->
            def cols = line.split('\\t')
            tuple(row_id, row_json, screen_json, n_loci, cols[0] as Integer, cols[1], cols[2])
        }
    }

    PREP_LOCUS(locus_tasks_ch)

    method_tasks_ch = PREP_LOCUS.out.flatMap { row_id, row_json, screen_json, n_loci, region_id, region_type, scratch_loci_root ->
        methods_list.collect { m_id ->
            tuple(row_id, row_json, screen_json, n_loci, region_id, region_type, scratch_loci_root, m_id)
        }
    }

    method_split_ch = method_tasks_ch.branch { row_id, row_json, screen_json, n_loci, region_id, region_type, scratch_loci_root, m_id ->
        light: m_id == 'p_ld'
            return tuple(row_id, row_json, screen_json, n_loci, region_id, region_type, scratch_loci_root)
        heavy: true
            return tuple(row_id, row_json, screen_json, n_loci, region_id, region_type, scratch_loci_root, m_id)
    }

    RUN_P_LD(method_split_ch.light)
    RUN_BAYESIAN_METHOD(method_split_ch.heavy)

    locus_done_ch = RUN_P_LD.out
        .mix(RUN_BAYESIAN_METHOD.out)
        .map { row_id, row_json, screen_json, n_loci, region_id, scratch_loci_root ->
            tuple(groupKey("${row_id}::${region_id}", methods_list.size()), row_id, row_json, screen_json, n_loci, region_id, scratch_loci_root)
        }
        .groupTuple(by: 0)
        .map { locus_gkey, row_ids, row_jsons, screen_jsons, n_locis, region_ids, loci_roots ->
            tuple(row_ids[0], row_jsons[0], screen_jsons[0], n_locis[0], region_ids[0], loci_roots[0])
        }

    CLEAN_LOCUS_SCRATCH(locus_done_ch)

    completed_loci_ch = CLEAN_LOCUS_SCRATCH.out
        .map { row_id, row_json, screen_json, n_loci, scratch_loci_root ->
            tuple(groupKey(row_id, n_loci as Integer), row_id, row_json, screen_json, scratch_loci_root)
        }
        .groupTuple(by: 0)
        .map { gwas_gkey, row_ids, row_jsons, screen_jsons, loci_roots ->
            tuple(row_ids[0], row_jsons[0], screen_jsons[0], loci_roots[0])
        }

    to_finalise_ch = branched_ch.zero_loci.mix(completed_loci_ch)
    FINALISE_GWAS(to_finalise_ch)

    UPDATE_STATE(FINALISE_GWAS.out.records.collect().ifEmpty([]))
}
