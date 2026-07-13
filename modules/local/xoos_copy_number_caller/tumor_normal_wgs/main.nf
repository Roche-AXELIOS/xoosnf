process XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/copy_number_caller:1.1.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(small_variant_vcf), path(small_variant_vcf_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(mappability_bigwig)
    tuple val(meta4), path(cnv_seed_file)
    tuple val(meta5), path(blocklist_bed)

    output:
    tuple val(meta), path("copy_number_caller/${prefix}_somatic_copy_number_callset.vcf.gz"), emit: vcf
    tuple val(meta), path("copy_number_caller/*"), emit: results_dir
    tuple val("${task.process}"), val('copy_number_caller'), eval("copy_number_caller --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    mkdir -p "copy_number_caller"

    copy_number_caller tumor-normal-wgs \\
        --tumor-bam '${tumor_bam}' \\
        --normal-bam '${normal_bam}' \\
        --tumor-sample-name '${meta.tumor_id}' \\
        --normal-sample-name '${meta.normal_id}' \\
        --vcf '${small_variant_vcf}' \\
        --reference '${reference}' \\
        --seed-segments '${cnv_seed_file ?: '/resources/hg38_seed.seg'}' \\
        --mappability-bigwig '${mappability_bigwig ?: '/resources/hg38_mappability.bw'}' \\
        --blocklist-bed '${blocklist_bed ?: '/resources/hg38_blocklist.bed'}' \\
        --output-dir 'copy_number_caller' \\
        ${args}

    for file in copy_number_caller/*; do
        mv -- "\$file" "copy_number_caller/${prefix}_\$(basename "\$file")"
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "copy_number_caller"
    touch 'copy_number_caller/${prefix}_somatic_copy_number_callset.vcf.gz'
    touch 'copy_number_caller/${prefix}_somatic_copy_number_callset.seg'
    touch 'copy_number_caller/${prefix}_purity_ploidy_grid.tsv'
    touch 'copy_number_caller/${prefix}_igv_visualization.xml'
    touch 'copy_number_caller/${prefix}_metrics.tsv'
    touch 'copy_number_caller/${prefix}_augmented_intervals.bed'
    touch 'copy_number_caller/${prefix}_normal_coverage.bed'
    touch 'copy_number_caller/${prefix}_normal_corrected_coverage.bed'
    touch 'copy_number_caller/${prefix}_tumor_coverage.bed'
    touch 'copy_number_caller/${prefix}_tumor_corrected_coverage.bed'
    touch 'copy_number_caller/${prefix}_log_ratios.bw'
    touch 'copy_number_caller/${prefix}_log_ratios.bed'
    touch 'copy_number_caller/${prefix}_log_ratio_segments.bed'
    touch 'copy_number_caller/${prefix}_mapping_quality.bed'
    touch 'copy_number_caller/${prefix}_bafs.bed'
    touch 'copy_number_caller/${prefix}_bafs.bw'
    """
}
