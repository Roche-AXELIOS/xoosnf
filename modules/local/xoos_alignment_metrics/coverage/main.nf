process XOOS_ALIGNMENT_METRICS_COVERAGE {
    conda false
    container 'ghcr.io/roche-axelios/xoos/alignment_metrics:1.1.0'

    label 'process_high'
    tag "${meta1.id}"

    input:
    tuple val(meta1), path(bam), path(bam_index)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta1), path("coverage"), optional: true, emit: coverage_metrics
    tuple val("${task.process}"), val('alignment_metrics'), eval('alignment_metrics --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reference_opt = fasta ? "--reference ${fasta}" : ""
    def bed_opt = bed ? "--bed-input ${bed}" : ""
    def args = task.ext.args ?: ''
    """
    alignment_metrics coverage \\
        --bam-input '${bam}' \\
        ${bed_opt} \\
        ${reference_opt} \\
        --output-dir . \\
        ${args}
    """

    stub:
    """
    install -D /dev/null 'coverage/coverage_stats.tsv'
    """
}
