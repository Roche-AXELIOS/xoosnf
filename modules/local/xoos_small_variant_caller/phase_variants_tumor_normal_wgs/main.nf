process XOOS_SMALL_VARIANT_CALLER_PHASE_VARIANTS_TUMOR_NORMAL_WGS {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.2.0'

    input:
    tuple val(meta), path(bam), path(bam_index), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("${prefix}_phased.vcf.gz"), path("${prefix}_phased.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('phase_variants'), eval("phase_variants --version"), topic: versions, emit: versions_phase_variants

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: meta.id
    """
    phase_variants tumor-normal-wgs \\
        --bam-input '${bam}' \\
        --vcf-input '${vcf}' \\
        --output-dir '.' \\
        --vcf-output '${prefix}_phased.vcf.gz' \\
        --tumor-sample-name '${meta.tumor_sample_name ?: meta.tumor_id ?: meta.id}' \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    gzip -n < /dev/null > '${prefix}_phased.vcf.gz'
    touch '${prefix}_phased.vcf.gz.tbi'
    """
}
