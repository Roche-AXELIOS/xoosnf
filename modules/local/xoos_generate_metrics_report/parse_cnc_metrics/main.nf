process XOOS_PARSE_CNC_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.1.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(results_files)

    output:
    tuple val(meta), path("${meta.id}.copy_number_caller.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    metrics_file=\$(find . -maxdepth 1 -name '*_metrics.tsv' | head -1)
    if [ -z "\$metrics_file" ]; then
        echo "ERROR: *_metrics.tsv not found in staged files" >&2
        exit 1
    fi

    generate-metrics-report parse-cnc-metrics \\
        --sample-name '${meta.id}' \\
        --input "\$metrics_file" \\
        --output '${meta.id}.copy_number_caller.json'
    """

    stub:
    """
    echo '[]' > '${meta.id}.copy_number_caller.json'
    """
}
