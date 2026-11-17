nextflow.enable.dsl = 2

// Mirrors GWASprocessing/pipeline/main.nf:
// 1. PREPARE_MANIFEST discovers/validates inputs, skips complete rows via .finemap_state.tsv, and emits rows/*.json
//    so phenotype labels and file paths are never interpolated into shell strings.
// 2. FINEMAP_ROW processes one GWAS per task with fault isolation and copies task logs to <outdir>/<study>/logs/.
// 3. UPDATE_STATE updates .finemap_state.tsv and renders <outdir>/manifest.html.

def opt(String flag, def val) {
    return (val != null && val.toString().trim() != '') ? "${flag} '${val}'" : ""
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
      ${opt('--local-ld-panel', params.local_ld_panel)} \\
      --cache-dir '${params.cache_dir}' \\
      --population '${params.population}' \\
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

process FINEMAP_ROW {
    tag "${row_json.baseName}"
    cache false
    afterScript "python3 '${projectDir}/scripts/prepare_manifest.py' copy-logs --outdir '${params.outdir}' --row-json '${row_json}' || true"

    input:
    path row_json

    output:
    path '*.row_record.json', emit: records

    script:
    """
    python3 '${projectDir}/scripts/run_finemapping.py' \\
      --row-json '${row_json}' \\
      --output-record '${row_json.baseName}.row_record.json'
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

    def chain_file = new File(params.chain.toString()).canonicalFile
    def fasta_file = new File(params.grch37_fasta.toString()).canonicalFile
    def fai_file = new File(params.grch37_fai.toString()).canonicalFile
    [chain_file, fasta_file, fai_file].each { f ->
        if (!f.isFile()) error "Reference file does not exist: ${f}"
    }
    if (new File("${fasta_file}.fai").canonicalFile != fai_file) {
        error "--grch37_fai must be the .fai index adjacent to --grch37_fasta"
    }

    def output_dir = new File(params.outdir.toString()).canonicalFile
    params.outdir = output_dir.toString()
    if (!output_dir.exists() && !output_dir.mkdirs()) {
        error "Could not create output directory: ${output_dir}"
    }

    PREPARE_MANIFEST()
    def row_records = FINEMAP_ROW(PREPARE_MANIFEST.out.rows.flatten())
    UPDATE_STATE(row_records.records.collect().ifEmpty([]))
}
