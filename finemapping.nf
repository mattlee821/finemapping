nextflow.enable.dsl = 2

process FINEMAP_ONE_GWAS {
    tag "${study_id}/${gwas_id}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    cpus 2
    memory '8 GB'
    time '48h'

    input:
    tuple val(study_id), val(gwas_id),
          path(vcf), path(vcf_index), path(metadata), path(qc),
          path(chain), path(grch37_fasta), path(grch37_fai),
          path(gene_annotation)

    output:
    path "${study_id}/${gwas_id}", emit: results

    script:
    def reference_size_arg = params.reference_sample_size ?
        "--reference-sample-size ${params.reference_sample_size}" : ""
    """
    python3 '${projectDir}/scripts/run_finemapping.py' \\
      --vcf '${vcf}' \\
      --metadata '${metadata}' \\
      --qc '${qc}' \\
      --chain '${chain}' \\
      --grch37-fasta '${grch37_fasta}' \\
      --gene-annotation '${gene_annotation}' \\
      --output-root . \\
      --cache-dir '${params.cache_dir}' \\
      --population '${params.population}' \\
      ${reference_size_arg} \\
      --significance-lp '${params.significance_lp}' \\
      --merge-distance-bp '${params.merge_distance_bp}' \\
      --region-flank-bp '${params.region_flank_bp}' \\
      --gene-window-bp '${params.gene_window_bp}' \\
      --max-loci '${params.max_loci}'
    """
}

workflow {
    def required = [
        'vcf', 'vcf_index', 'metadata', 'qc', 'chain', 'grch37_fasta',
        'grch37_fai', 'gene_annotation', 'study_id', 'gwas_id', 'outdir'
    ]
    def missing = required.findAll { !params[it] }
    if (missing) {
        error "Missing required parameter(s): ${missing.collect { "--${it}" }.join(', ')}"
    }

    def input_paths = [
        params.vcf, params.vcf_index, params.metadata, params.qc, params.chain,
        params.grch37_fasta, params.grch37_fai, params.gene_annotation
    ].collect { new File(it).canonicalFile }
    def absent = input_paths.findAll { !it.isFile() }
    if (absent) {
        error "Input file(s) do not exist: ${absent.join(', ')}"
    }
    def expected_fai = new File("${params.grch37_fasta}.fai").canonicalFile
    if (!expected_fai.isFile() || expected_fai != input_paths[6]) {
        error "--grch37_fai must be the index for --grch37_fasta"
    }
    if (params.max_loci.toString().toInteger() < 0) {
        error "--max_loci must be non-negative"
    }
    def output_dir = new File(params.outdir).canonicalFile
    if (input_paths.any { it.toPath().startsWith(output_dir.toPath()) }) {
        error "Output directory must not contain any input file: ${output_dir}"
    }
    params.outdir = output_dir.toString()
    if (!output_dir.exists() && !output_dir.mkdirs()) {
        error "Could not create output directory: ${output_dir}"
    }

    def inputs = Channel.of(tuple(
        params.study_id.toString(),
        params.gwas_id.toString(),
        file(params.vcf),
        file(params.vcf_index),
        file(params.metadata),
        file(params.qc),
        file(params.chain),
        file(params.grch37_fasta),
        file(params.grch37_fai),
        file(params.gene_annotation)
    ))
    FINEMAP_ONE_GWAS(inputs)
}
