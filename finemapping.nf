nextflow.enable.dsl = 2

process FINEMAP_ONE_GWAS {
    tag "${study_id}/${gwas_id}"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(study_id), val(gwas_id), val(trait_label), val(ancestry_label), val(gene_id_label),
          path(vcf), path(vcf_index),
          val(has_metadata), path(metadata),
          val(has_qc), path(qc),
          path(chain), path(grch37_fasta), path(grch37_fai),
          val(has_gene_annot), path(gene_annotation)

    output:
    path "${study_id}/${gwas_id}", emit: results

    script:
    def meta_arg = has_metadata ? "--metadata '${metadata}'" : ""
    def qc_arg = has_qc ? "--qc '${qc}'" : ""
    def gene_annot_arg = has_gene_annot ? "--gene-annotation '${gene_annotation}'" : ""
    def trait_arg = trait_label ? "--trait '${trait_label}'" : ""
    def anc_arg = ancestry_label ? "--ancestry '${ancestry_label}'" : ""
    def gid_arg = gene_id_label ? "--gene-id '${gene_id_label}'" : ""
    def panukb_arg = params.panukb_dir ? "--panukb-dir '${params.panukb_dir}'" : ""
    def local_ld_arg = params.local_ld_panel ? "--local-ld-panel '${params.local_ld_panel}'" : ""
    def ref_n_arg = params.reference_sample_size ? "--reference-sample-size ${params.reference_sample_size}" : ""
    def eff_n_arg = params.effective_n ? "--effective-n ${params.effective_n}" : ""
    def cis_lp_arg = params.cis_significance_lp ? "--cis-significance-lp ${params.cis_significance_lp}" : ""
    def mhc_arg = params.exclude_mhc.toString().toBoolean() ? "--exclude-mhc" : ""
    def keep_arg = params.keep_intermediates.toString().toBoolean() ? "--keep-intermediates" : ""
    """
    python3 '${projectDir}/scripts/run_finemapping.py' \\
      --vcf '${vcf}' \\
      --study-id '${study_id}' \\
      --gwas-id '${gwas_id}' \\
      ${meta_arg} \\
      ${qc_arg} \\
      ${trait_arg} \\
      ${anc_arg} \\
      ${gid_arg} \\
      --chain '${chain}' \\
      --grch37-fasta '${grch37_fasta}' \\
      ${gene_annot_arg} \\
      ${panukb_arg} \\
      ${local_ld_arg} \\
      --output-root . \\
      --cache-dir '${params.cache_dir}' \\
      --population '${params.population}' \\
      ${ref_n_arg} \\
      ${eff_n_arg} \\
      --significance-lp '${params.significance_lp}' \\
      ${cis_lp_arg} \\
      --region-flank-bp '${params.region_flank_bp}' \\
      --max-region-width-bp '${params.max_region_width_bp}' \\
      --merge-distance-bp '${params.merge_distance_bp}' \\
      --gene-window-bp '${params.gene_window_bp}' \\
      ${mhc_arg} \\
      ${keep_arg} \\
      --max-loci '${params.max_loci}' \\
      --max-causal '${params.max_causal}' \\
      --max-iter '${params.max_iter}'
    """
}

process RENDER_MANIFEST_REPORT {
    tag "manifest_report"

    input:
    val(completed_dirs)
    val(outdir_str)

    output:
    val(true), emit: done

    script:
    """
    python3 '${projectDir}/scripts/render_manifest_report.py' --outdir '${outdir_str}'
    """
}

workflow {
    def ref_required = ['chain', 'grch37_fasta', 'grch37_fai', 'outdir']
    def missing_ref = ref_required.findAll { !params[it] }
    if (missing_ref) {
        error "Missing required parameter(s): ${missing_ref.collect { "--${it}" }.join(', ')}"
    }
    if (!params.manifest && !params.vcf) {
        error "Supply either --vcf <path> (with --study_id and --gwas_id) or --manifest <tsv>"
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

    def has_gene = params.gene_annotation != null && params.gene_annotation.toString().trim() != ''
    def gene_file = has_gene ? new File(params.gene_annotation.toString()).canonicalFile : chain_file
    if (has_gene && !gene_file.isFile()) {
        error "--gene_annotation file does not exist: ${gene_file}"
    }

    def output_dir = new File(params.outdir.toString()).canonicalFile
    params.outdir = output_dir.toString()
    if (!output_dir.exists() && !output_dir.mkdirs()) {
        error "Could not create output directory: ${output_dir}"
    }

    def records = []
    if (params.manifest) {
        def mf = new File(params.manifest.toString()).canonicalFile
        if (!mf.isFile()) error "Manifest file does not exist: ${mf}"
        def lines = mf.readLines().findAll { it.trim() && !it.trim().startsWith('#') }
        if (lines.size() < 2) error "Manifest TSV has no data rows: ${mf}"
        def header = lines[0].split('\t', -1).collect { it.trim() }
        lines.drop(1).eachWithIndex { line, idx ->
            def vals = line.split('\t', -1)
            def row = [:]
            header.eachWithIndex { col, c_i -> row[col] = c_i < vals.size() ? vals[c_i].trim() : '' }
            def vcf_p = row.vcf ?: row.vcf_path
            if (!vcf_p) error "Manifest row ${idx + 2} lacks 'vcf' column"
            def vf = new File(vcf_p).canonicalFile
            if (!vf.isFile()) error "VCF in manifest row ${idx + 2} not found: ${vf}"
            def idx_p = row.vcf_index ?: (new File("${vf}.tbi").isFile() ? "${vf}.tbi" : "${vf}.csi")
            def vidx = new File(idx_p).canonicalFile
            if (!vidx.isFile()) error "Tabix/CSI index missing for ${vf}"

            def s_id = row.study_id ?: params.study_id ?: 'study'
            def g_id = row.gwas_id ?: row.output_id ?: vf.name.replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
            def meta_cand = row.metadata ? new File(row.metadata) : new File(vf.parentFile, 'metadata.json')
            def qc_cand = row.qc ? new File(row.qc) : new File(vf.parentFile, 'qc.json')
            def has_m = meta_cand.isFile()
            def has_q = qc_cand.isFile()

            records << tuple(
                s_id.toString(),
                g_id.toString(),
                (row.trait ?: '').toString(),
                (row.ancestry ?: '').toString(),
                (row.gene_id ?: row.uniprot_id ?: '').toString(),
                file(vf),
                file(vidx),
                has_m,
                file(has_m ? meta_cand.canonicalFile : vf),
                has_q,
                file(has_q ? qc_cand.canonicalFile : vf),
                file(chain_file),
                file(fasta_file),
                file(fai_file),
                has_gene,
                file(gene_file)
            )
        }
    } else {
        if (!params.study_id || !params.gwas_id) {
            error "Single-GWAS mode requires --study_id and --gwas_id"
        }
        def vf = new File(params.vcf.toString()).canonicalFile
        if (!vf.isFile()) error "--vcf does not exist: ${vf}"
        def idx_p = params.vcf_index ?: (new File("${vf}.tbi").isFile() ? "${vf}.tbi" : "${vf}.csi")
        def vidx = new File(idx_p.toString()).canonicalFile
        if (!vidx.isFile()) error "VCF index does not exist: ${vidx}"
        if (vf.toPath().startsWith(output_dir.toPath())) {
            error "Output directory must not contain input VCF: ${output_dir}"
        }
        def has_m = params.metadata != null && params.metadata.toString().trim() != ''
        def mf = has_m ? new File(params.metadata.toString()).canonicalFile : vf
        if (has_m && !mf.isFile()) error "--metadata does not exist: ${mf}"
        def has_q = params.qc != null && params.qc.toString().trim() != ''
        def qf = has_q ? new File(params.qc.toString()).canonicalFile : vf
        if (has_q && !qf.isFile()) error "--qc does not exist: ${qf}"

        records << tuple(
            params.study_id.toString(),
            params.gwas_id.toString(),
            (params.trait ?: '').toString(),
            (params.ancestry ?: '').toString(),
            (params.gene_id ?: '').toString(),
            file(vf),
            file(vidx),
            has_m,
            file(mf),
            has_q,
            file(qf),
            file(chain_file),
            file(fasta_file),
            file(fai_file),
            has_gene,
            file(gene_file)
        )
    }

    def out_ch = FINEMAP_ONE_GWAS(Channel.fromList(records))
    RENDER_MANIFEST_REPORT(out_ch.results.collect(), params.outdir)
}
