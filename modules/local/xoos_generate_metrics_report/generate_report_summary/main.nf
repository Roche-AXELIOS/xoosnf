def expandMultipleOptions(optionName, values) {
    return values.collect { value -> "${optionName} '${value}'" }.join(' \\\n')
}

process XOOS_GENERATE_REPORT_SUMMARY {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(metrics)
    val(summary_patterns)

    output:
    tuple val(meta), path('metrics-report-summary.tsv'), emit: summary_report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report generate-report \\
        ${expandMultipleOptions('--input', metrics)} \\
        ${expandMultipleOptions('--pattern', summary_patterns)} \\
        --output metrics-report-summary.tsv
    """

    stub:
    """
    echo -e "metric\tvalue" > metrics-report-summary.tsv
    """
}
