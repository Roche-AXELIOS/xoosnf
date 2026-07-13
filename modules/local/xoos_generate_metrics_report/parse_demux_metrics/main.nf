process XOOS_PARSE_DEMUX_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.1.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(run_stats), path(sample_stats)

    output:
    tuple val(meta), path('demux_metrics.json'), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-demux-metrics \\
        --run-stats '${run_stats}' \\
        --sample-stats '${sample_stats}' \\
        --output demux_metrics.json
    """

    stub:
    """
    echo '{}' > demux_metrics.json
    """
}
