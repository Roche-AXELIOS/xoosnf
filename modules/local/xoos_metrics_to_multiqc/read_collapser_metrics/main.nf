process XOOS_READ_COLLAPSER_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.2.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(read_collapser_metrics)

    output:
    tuple val(meta), path("*_read_collapser_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('xoos_to_multiqc.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p '${prefix}_read_collapser_multiqc_input'
    xoos_to_multiqc.py \\
        --sample-name ${prefix} \\
        --read-collapser-summary ${read_collapser_metrics}/metrics/summary_stats.tsv \\
        --read-collapser-cluster-sizes ${read_collapser_metrics}/metrics/cluster_size_distributions.tsv \\
        --output-dir '${prefix}_read_collapser_multiqc_input' \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p '${prefix}_read_collapser_multiqc_input'
    touch '${prefix}_read_collapser_multiqc_input/xoos_read_collapser_metrics_mqc.json'
    touch '${prefix}_read_collapser_multiqc_input/xoos_read_collapser_cluster_sizes_mqc.json'
    """
}
