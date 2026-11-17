nextflow.enable.dsl = 2

process BUILD_MANIFEST {
    tag 'manifest'
    cache false
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    val input_dir

    output:
    path 'manifest.tsv', emit: manifest
    path 'discovery_summary.json', emit: discovery

    script:
    """
    python3 '${projectDir}/scripts/gwas_pilot.py' manifest \\
      --input-dir '${input_dir}' \\
      --max-inputs ${params.max_inputs} \\
      --output manifest.tsv \\
      --summary discovery_summary.json \\
      --existing-manifest '${params.outdir}/manifest.tsv' \\
      --existing-summary '${params.outdir}/discovery_summary.json'
    """
}

process PREFLIGHT {
    tag "$study_id"
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true
    cpus 1
    memory '1 GB'
    time '10m'

    input:
    tuple val(study_id),
          path(vcf), val(source_vcf),
          path(index), val(source_index),
          path(metadata), val(source_metadata),
          path(qc), val(source_qc)

    output:
    path "${study_id}.preflight.json", emit: reports

    script:
    """
    python3 '${projectDir}/scripts/gwas_pilot.py' preflight \\
      --study-id '${study_id}' \\
      --vcf '${vcf}' \\
      --index '${index}' \\
      --metadata '${metadata}' \\
      --qc '${qc}' \\
      --source-vcf '${source_vcf}' \\
      --source-index '${source_index}' \\
      --source-metadata '${source_metadata}' \\
      --source-qc '${source_qc}' \\
      --output '${study_id}.preflight.json'
    """
}

process SUMMARIZE {
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path reports

    output:
    path 'pilot_summary.tsv'

    script:
    def report_args = reports.collect { "'${it}'" }.join(' ')
    """
    python3 '${projectDir}/scripts/gwas_pilot.py' summarize \\
      --reports ${report_args} \\
      --output pilot_summary.tsv
    """
}

workflow {
    if (!params.input_dir) {
        error "Missing required parameter: --input_dir"
    }
    if (!params.outdir) {
        error "Missing required parameter: --outdir"
    }
    def max_inputs
    def max_forks
    try {
        max_inputs = params.max_inputs.toString().toInteger()
        max_forks = params.max_forks.toString().toInteger()
    } catch (NumberFormatException ignored) {
        error "--max_inputs and --max_forks must be integers"
    }
    if (max_inputs < 0 || max_forks < 1) {
        error "--max_inputs must be non-negative and --max_forks must be positive"
    }

    def input_root = new File(params.input_dir).canonicalFile
    if (!input_root.isDirectory()) {
        error "Input directory does not exist: ${input_root}"
    }
    def output_dir = new File(params.outdir).canonicalFile
    params.outdir = output_dir.toString()
    if (output_dir.toPath().startsWith(input_root.toPath())) {
        error "Output directory must not be inside the input directory: ${output_dir}"
    }
    if (output_dir.exists() && !output_dir.isDirectory()) {
        error "Output path is not a directory: ${output_dir}"
    }
    def existing_count = 0
    if (output_dir.exists()) {
        def existing_outputs = output_dir.list()
        if (existing_outputs == null) {
            error "Could not read output directory: ${output_dir}"
        }
        existing_count = existing_outputs.length
    }
    if (existing_count > 0 && !workflow.resume) {
        error "Output directory must be empty or new: ${output_dir}"
    }
    if (existing_count > 0 && !new File(output_dir, 'manifest.tsv').isFile()) {
        error "Cannot resume without the original manifest in: ${output_dir}"
    }
    if (!output_dir.exists() && !output_dir.mkdirs()) {
        error "Could not create output directory: ${output_dir}"
    }

    BUILD_MANIFEST(Channel.value(input_root.toString()))
    def input_rows = BUILD_MANIFEST.out.manifest
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.study_id,
                file(row.vcf),
                row.vcf,
                file(row.index),
                row.index,
                file(row.metadata),
                row.metadata,
                file(row.qc),
                row.qc
            )
        }

    PREFLIGHT(input_rows)
    SUMMARIZE(PREFLIGHT.out.reports.collect())
}
