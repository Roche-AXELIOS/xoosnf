process XOOS_READ_COLLAPSER_RESCUE_SECONDARY {
    conda false
    container 'ghcr.io/roche-axelios/xoos/read_collapser:1.1.0'

    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: rescued_bam
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
        --prefix ${meta.id} \\
        ${args}
    mv rescue_secondary_metrics/${meta.id}.rescued_secondary.bam ${meta.id}.bam
    """

    stub:
    """
    mkdir -p 'rescue_secondary_metrics'
    touch '${meta.id}.bam'
    touch 'rescue_secondary_metrics/rescue_secondary_metrics.tsv'
    """
}
