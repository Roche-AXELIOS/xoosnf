process XOOS_COPY_NUMBER_CALLER_GERMLINE_WGS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/copy_number_caller:1.1.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index), path(small_variant_vcf), path(small_variant_vcf_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(mappability_bigwig)
    tuple val(meta4), path(cnv_seed_file)

    output:
    tuple val(meta), path("copy_number_caller/${prefix}_germline_copy_number_callset.vcf.gz"), emit: vcf
    tuple val(meta), path("copy_number_caller/*"), emit: results_dir
    tuple val("${task.process}"), val('copy_number_caller'), eval("copy_number_caller --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def vcf_arg = small_variant_vcf ? "--vcf '${small_variant_vcf}'" : ''
    """
    mkdir -p "copy_number_caller"

    copy_number_caller germline-wgs \\
        --bam '${bam}' \\
        --reference '${reference}' \\
        --mappability-bigwig '${mappability_bigwig ?: '/resources/hg38_mappability.bw'}' \\
        --seed-segments '${cnv_seed_file ?: '/resources/hg38_seed.seg'}' \\
        --sample-name '${meta.id}' \\
        --output-dir 'copy_number_caller' \\
        ${vcf_arg} \\
        ${args}

    for file in copy_number_caller/*; do
        mv -- "\$file" "copy_number_caller/${prefix}_\$(basename "\$file")"
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "copy_number_caller"
    touch 'copy_number_caller/${prefix}_germline_copy_number_callset.seg'
    gzip -n < /dev/null > 'copy_number_caller/${prefix}_germline_copy_number_callset.vcf.gz'
    touch 'copy_number_caller/${prefix}_germline_copy_number_callset.vcf.gz.tbi'
    touch 'copy_number_caller/${prefix}_igv_visualization.xml'
    touch 'copy_number_caller/${prefix}_metrics.tsv'
    touch 'copy_number_caller/${prefix}_augmented_intervals.bed'
    touch 'copy_number_caller/${prefix}_coverage.bed'
    touch 'copy_number_caller/${prefix}_corrected_coverage.bed'
    touch 'copy_number_caller/${prefix}_log_ratios.bw'
    touch 'copy_number_caller/${prefix}_log_ratios.bed'
    touch 'copy_number_caller/${prefix}_mapping_quality.bed'
    touch 'copy_number_caller/${prefix}_bafs.bed'
    touch 'copy_number_caller/${prefix}_bafs.bw'
    """
}
