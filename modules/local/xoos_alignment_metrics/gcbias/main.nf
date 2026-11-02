process XOOS_ALIGNMENT_METRICS_GC_BIAS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/alignment_metrics:1.2.0'

    label 'process_medium'
    tag "${meta1.id}"

    input:
    tuple val(meta1), path(bam), path(bam_index)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta1), path('gc_bias'), emit: gc_bias_metrics
    tuple val("${task.process}"), val('alignment_metrics'), eval('alignment_metrics --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    alignment_metrics gcbias \\
        --bam '${bam}' \\
        --reference '${fasta}' \\
        --output-dir . \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    install -D /dev/null 'gc_bias/gc_bias.tsv'
    install -D /dev/null 'gc_bias/gc_bias_summary.tsv'
    """
}
