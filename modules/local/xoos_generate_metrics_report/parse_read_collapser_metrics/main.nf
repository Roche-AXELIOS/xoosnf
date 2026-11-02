process XOOS_PARSE_READ_COLLAPSER_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(read_collapser_metrics)

    output:
    tuple val(meta), path("${meta.id}.read_collapser.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-read-collapser-metrics \\
        --sample-name '${meta.id}' \\
        --source 'read_collapser' \\
        --input '${read_collapser_metrics}/metrics/summary_stats.tsv' \\
        --input '${read_collapser_metrics}/metrics/cluster_size_distribution_summary.tsv' \\
        --output '${meta.id}.read_collapser.json'
    """

    stub:
    """
    echo '{}' > '${meta.id}.read_collapser.json'
    """
}
