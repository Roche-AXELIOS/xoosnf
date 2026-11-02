process XOOS_PARSE_GC_BIAS_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'
    conda false

    label 'process_single'
    tag "${meta.id}, ${source}"

    input:
    tuple val(meta), val(source), path(gc_bias)

    output:
    tuple val(meta), path("${meta.id}.gc_bias_${source}.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate-metrics-report parse-gc-bias \\
        --sample-name '${meta.id}' \\
        --source 'gc_bias/${source}' \\
        --summary '${gc_bias}/${prefix}.gc_bias_summary.tsv' \\
        --detail '${gc_bias}/${prefix}.gc_bias.tsv' \\
        --output '${meta.id}.gc_bias_${source}.json'
    """

    stub:
    """
    echo '{}' > '${meta.id}.gc_bias_${source}.json'
    """
}
