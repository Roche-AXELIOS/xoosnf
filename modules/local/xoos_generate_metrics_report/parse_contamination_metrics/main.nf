process XOOS_PARSE_CONTAMINATION_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(results_dir)

    output:
    tuple val(meta), path("${meta.id}.contamination_detection.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-contamination-metrics \\
        --sample-name '${meta.id}' \\
        --input '${results_dir}/contamination_estimation_results.tsv' \\
        --output '${meta.id}.contamination_detection.json'
    """

    stub:
    """
    echo '[]' > '${meta.id}.contamination_detection.json'
    """
}
