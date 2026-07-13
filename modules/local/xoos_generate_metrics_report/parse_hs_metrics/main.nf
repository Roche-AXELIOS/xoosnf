process XOOS_PARSE_HS_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.1.0'
    conda false

    label 'process_single'
    tag "${meta.id}, ${source}"

    input:
    tuple val(meta), val(source), path(hs_metrics)

    output:
    tuple val(meta), path("${meta.id}.hs_metrics_${source}.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-hs-metrics \\
        --sample-name '${meta.id}' \\
        --source 'hs_metrics/${source}' \\
        --input '${hs_metrics}' \\
        --output '${meta.id}.hs_metrics_${source}.json'
    """

    stub:
    """
    echo '{}' > '${meta.id}.hs_metrics_${source}.json'
    """
}
