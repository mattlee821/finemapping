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
      --min-maf '${params.min_maf}' \\
      --region-flank-bp '${params.region_flank_bp}' \\
      --max-region-width-bp '${params.max_region_width_bp}' \\
      --gene-window-bp '${params.gene_window_bp}' \\
      --exclude-mhc '${params.exclude_mhc}' \\
      --max-loci '${params.max_loci}' \\
      --max-causal '${params.max_causal}' \\
      --max-iter '${params.max_iter}' \\
      --cis-max-iter '${params.cis_max_iter}' \\
      --susie-refine '${params.susie_refine}' \\
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
    afterScript { "python3 '${projectDir}/scripts/prepare_manifest.py' copy-logs --outdir '${params.outdir}' --row-json '${row_json}' --step-name '01_screen_gwas' || true" }

    input:
    tuple val(row_id), path(row_json)

    output:
    tuple val(row_id), path(row_json), path('screen_summary.json'), path('loci_list.tsv')

    script:
    def tmp_export = params.temp_dir ? "mkdir -p '${params.temp_dir}' && export TMPDIR='${params.temp_dir}'" : ""
    """
    ${tmp_export}
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
    "--min-maf", str(row.get("min_maf", 0.01)),
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
    w.writerow(["locus_index", "region_id", "region_type", "mem_gb"])
    for idx, reg in enumerate(summary.get("selected_regions") or []):
        w.writerow([idx, reg["region_id"], reg["region_type"], reg.get("recommended_mem_gb", 16)])
PY
    """
}

process FINEMAP_LOCUS {
    tag "${row_id}:${region_id}"
    label 'method_heavy'
    cache false
    memory {
        def base = (mem_gb ?: 16) as Integer
        def mult = task.attempt == 1 ? 1.0 : (task.attempt == 2 ? 1.75 : 2.5)
        def req = Math.min(350, Math.max(base, Math.ceil(base * mult) as Integer))
        "${req} GB"
    }
    afterScript { "if [ -s .command.err ]; then python3 '${projectDir}/scripts/prepare_manifest.py' copy-logs --outdir '${params.outdir}' --row-json '${row_json}' --step-name '02_locus_worker_${region_id}' || true; fi" }

    input:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), val(locus_index), val(region_id), val(region_type), val(mem_gb)

    output:
    tuple val(row_id), path(row_json), path(screen_json), val(n_loci), env('SCRATCH_LOCI_ROOT')

    script:
    def scratch_root = params.scratch_dir ?: "${params.cache_dir}/scratch"
    """
    python3 '${projectDir}/scripts/locus_worker.py' \\
      --project-dir '${projectDir}' \\
      --row-json '${row_json}' \\
      --screen-json '${screen_json}' \\
      --locus-index '${locus_index}' \\
      --region-id '${region_id}' \\
      --region-type '${region_type}' \\
      --scratch-root '${scratch_root}' \\
      --temp-dir '${params.temp_dir ?: ""}' \\
      --bm-disk-cache '${params.bm_cache_dir ?: ""}' \\
      --attempt '${task.attempt}' \\
      --max-attempts '${(task.maxRetries ?: 0) + 1}' \\
      --out-root-file .scratch_loci_root
    SCRATCH_LOCI_ROOT="\$(cat .scratch_loci_root)"
    """
}

process FINALISE_GWAS {
    tag "${row_id}"
    cache false
    afterScript { "python3 '${projectDir}/scripts/prepare_manifest.py' copy-logs --outdir '${params.outdir}' --row-json '${row_json}' --step-name '04_finalise_gwas' || true" }

    input:
    tuple val(row_id), path(row_json), path(screen_json), val(scratch_loci_root)

    output:
    path '*.row_record.json', emit: records

    script:
    def tmp_export = params.temp_dir ? "mkdir -p '${params.temp_dir}' && export TMPDIR='${params.temp_dir}'" : ""
    """
    ${tmp_export}
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
    "--max-iter", str(row.get("max_iter", 100)),
    "--cis-max-iter", str(row.get("cis_max_iter", 200)),
    "--min-maf", str(row.get("min_maf", 0.01)),
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
# Rows reaching FINALISE_GWAS have already passed PREPARE_MANIFEST idempotency checks,
# so overwrite any partial non-log files left by a previously interrupted run.
cmd.append("--force")

try:
    subprocess.run(cmd, check=True)
    subprocess.run([
        sys.executable,
        str(Path(fin_script).with_name("prepare_manifest.py")),
        "update-state",
        "--state-file", str(Path(row["outdir"]) / ".finemap_state.tsv"),
        "--records", f"{row_id}.row_record.json",
        "--skip-manifest-report",
    ], check=False)
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
            def mem_gb = cols.size() > 3 ? (cols[3] as Integer) : 16
            tuple(row_id, row_json, screen_json, n_loci, cols[0] as Integer, cols[1], cols[2], mem_gb)
        }
    }

    FINEMAP_LOCUS(locus_tasks_ch)

    completed_loci_ch = FINEMAP_LOCUS.out
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
