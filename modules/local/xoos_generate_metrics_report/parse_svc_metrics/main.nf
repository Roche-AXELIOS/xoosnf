process XOOS_PARSE_SVC_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'
    conda false

    label 'process_single'
    tag "${meta.id}, ${source}"

    input:
    tuple val(meta), path(metrics_tsv), val(source)

    output:
    tuple val(meta), path("${meta.id}.${source}.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-svc-metrics \\
        --sample-name '${meta.id}' \\
        --source '${source}' \\
        --input '${metrics_tsv}' \\
        --output '${meta.id}.${source}.json'
    """

    stub:
    """
    echo '[]' > '${meta.id}.${source}.json'
    """
}
