process PICARD_COLLECT_GC_BIAS_METRICS {
    conda false
    container "ghcr.io/roche-csi/picard:3.4.0-roche.1"

    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta2), path(reference), path(reference_fai)

    output:
    tuple val(meta), path("CollectGcBiasMetrics"), emit: gc_bias
    tuple val("${task.process}"), val('picard'), val("3.4.0"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p 'CollectGcBiasMetrics'
    picard ${args2} CollectGcBiasMetrics \\
        -I ${bam} \\
        -O 'CollectGcBiasMetrics/${prefix}.gc_bias.tsv' \\
        -S 'CollectGcBiasMetrics/${prefix}.gc_bias_summary.tsv' \\
        -CHART 'CollectGcBiasMetrics/${prefix}.gc_bias.pdf' \\
        -R ${reference} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p CollectGcBiasMetrics
    touch CollectGcBiasMetrics/${prefix}.gc_bias.tsv
    touch CollectGcBiasMetrics/${prefix}.gc_bias_summary.tsv
    touch CollectGcBiasMetrics/${prefix}.gc_bias.pdf
    """
}
