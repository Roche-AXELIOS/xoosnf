def expandMultipleOptions(optionName, values) {
    return values.collect { value -> "${optionName} '${value}'" }.join(' \\\n')
}

process XOOS_GENERATE_REPORT {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.1.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(metrics)

    output:
    tuple val(meta), path('metrics-report.tsv'), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report generate-report \\
        ${expandMultipleOptions('--input', metrics)} \\
        --output metrics-report.tsv
    """

    stub:
    """
    echo -e "metric\tvalue" > metrics-report.tsv
    """
}
