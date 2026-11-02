process XOOS_COPY_NUMBER_CALLER_TUMOR_ONLY_WGS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/copy_number_caller:1.2.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(small_variant_vcf), path(small_variant_vcf_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(mappability_bigwig)
    tuple val(meta4), path(cnv_seed_file)
    tuple val(meta5), path(blocklist_bed)
    tuple val(meta6), path(population_vcf), path(population_vcf_index)
    tuple val(meta7), path(snp_blocklist_bed)

    output:
    tuple val(meta), path("copy_number_caller/${prefix}_somatic_copy_number_callset.vcf.gz"), emit: vcf
    tuple val(meta), path("copy_number_caller/${prefix}_somatic_copy_number_callset.seg"), emit: seg
    tuple val(meta), path("copy_number_caller/supplemental/${prefix}_purity_ploidy_grid.tsv"), emit: purity_ploidy
    tuple val(meta), path("copy_number_caller/*"), emit: results_dir
    tuple val("${task.process}"), val('copy_number_caller'), eval("copy_number_caller --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def mappability_arg = mappability_bigwig ? "--mappability-bigwig '${mappability_bigwig}'" : ''
    """
    mkdir -p "copy_number_caller"

    copy_number_caller tumor-only-wgs \\
        --tumor-bam '${tumor_bam}' \\
        --tumor-sample-name '${meta.tumor_sample_name ?: meta.tumor_id ?: meta.id}' \\
        --vcf '${small_variant_vcf}' \\
        --population-vcf '${population_vcf}' \\
        --reference '${reference}' \\
        --seed-segments '${cnv_seed_file ?: '/resources/hg38_seed.seg'}' \\
        ${mappability_arg} \\
        --blocklist-bed '${blocklist_bed ?: '/resources/hg38_tumor_only_wgs_blocklist.bed'}' \\
        ${snp_blocklist_bed ? "--snp-blocklist-bed '${snp_blocklist_bed}'" : ''} \\
        --output-dir 'copy_number_caller' \\
        --output-prefix '${prefix}' \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def metrics_dir = 'copy_number_caller/metrics'
    def tracks_dir = 'copy_number_caller/tracks'
    def supplemental_dir = 'copy_number_caller/supplemental'
    """
    mkdir -p '${metrics_dir}' '${tracks_dir}' '${supplemental_dir}'
    gzip -n < /dev/null > 'copy_number_caller/${prefix}_somatic_copy_number_callset.vcf.gz'
    touch 'copy_number_caller/${prefix}_somatic_copy_number_callset.vcf.gz.tbi'
    touch 'copy_number_caller/${prefix}_somatic_copy_number_callset.seg'
    touch '${supplemental_dir}/${prefix}_purity_ploidy_grid.tsv'
    touch '${metrics_dir}/${prefix}_metrics.tsv'
    touch '${tracks_dir}/${prefix}_igv_visualization.xml'
    touch '${tracks_dir}/${prefix}_log_ratios.bw'
    touch '${tracks_dir}/${prefix}_bafs.bw'
    touch '${supplemental_dir}/${prefix}_augmented_intervals.bed'
    touch '${supplemental_dir}/${prefix}_tumor_coverage.bed'
    touch '${supplemental_dir}/${prefix}_tumor_corrected_coverage.bed'
    touch '${supplemental_dir}/${prefix}_log_ratios.bed'
    touch '${supplemental_dir}/${prefix}_log_ratio_segments.seg'
    touch '${supplemental_dir}/${prefix}_mapping_quality.bed'
    touch '${supplemental_dir}/${prefix}_bafs.bed'
    touch '${supplemental_dir}/${prefix}_segment_likelihood.diploid_bias.tsv'
    touch '${supplemental_dir}/${prefix}_segment_likelihood.uniform.tsv'
    """
}
