process XOOS_ALIGNMENT_METRICS_ALL {
    conda false
    container 'ghcr.io/roche-axelios/xoos/alignment_metrics:1.1.0'

    label 'process_high'
    tag "${meta1.id}"

    input:
    tuple val(meta1), path(bam), path(bam_index)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta1), path("coverage"), emit: coverage_metrics
    tuple val(meta1), path("accuracy"), emit: accuracy_metrics
    tuple val(meta1), path("read_metrics"), emit: read_metrics
    tuple val("${task.process}"), val('alignment_metrics'), eval('alignment_metrics --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reference_opt = fasta ? "--reference ${fasta}" : ""
    def bed_opt = bed ? "--bed-input ${bed}" : ""
    def args = task.ext.args ?: ''
    """
    alignment_metrics all \\
        --bam-input '${bam}' \\
        ${bed_opt} \\
        ${reference_opt} \\
        --output-dir . \\
        ${args}
    """

    stub:
    """
    install -D /dev/null 'coverage/coverage_stats.tsv'
    install -D /dev/null 'accuracy/base_level_accuracy_summary.tsv'
    install -D /dev/null 'read_metrics/read_metrics_summary.tsv'
    """
}
