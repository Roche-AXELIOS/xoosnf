process XOOS_READ_COLLAPSER_RESCUE_SECONDARY {
    conda false
    container 'ghcr.io/roche-axelios/xoos/read_collapser:1.2.0'

    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index)

    output:
    tuple val(meta), path("${meta.id}_rescued.bam"), emit: rescued_bam
    tuple val(meta), path("rescue_secondary_metrics"), emit: metrics_dir
    tuple val("${task.process}"), val('read_collapser'), eval("read_collapser --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    read_collapser_rescue_secondary \\
        --bam ${bam} \\
        --output-dir rescue_secondary_metrics \\
        --output-prefix ${meta.id} \\
        --threads ${task.cpus} \\
        ${args}
    # The default output BAM is '<prefix>_rescued_secondary.bam' and metrics are
    # written under a 'metrics/' subdirectory as
    # '<prefix>_rescue_secondary_metrics.tsv'.
    mv rescue_secondary_metrics/${meta.id}_rescued_secondary.bam ${meta.id}_rescued.bam
    """

    stub:
    """
    mkdir -p 'rescue_secondary_metrics/metrics'
    touch '${meta.id}_rescued.bam'
    touch 'rescue_secondary_metrics/metrics/${meta.id}_rescue_secondary_metrics.tsv'
    """
}
