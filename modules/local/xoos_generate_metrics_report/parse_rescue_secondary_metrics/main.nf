process XOOS_PARSE_RESCUE_SECONDARY_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(rescue_secondary_metrics)

    output:
    tuple val(meta), path("${meta.id}.rescue_secondary.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-read-collapser-metrics \\
        --sample-name '${meta.id}' \\
        --source 'rescue_secondary' \\
        --input '${rescue_secondary_metrics}/metrics/${meta.id}_rescue_secondary_metrics.tsv' \\
        --output '${meta.id}.rescue_secondary.json'
    """

    stub:
    """
    echo '{}' > '${meta.id}.rescue_secondary.json'
    """
}
