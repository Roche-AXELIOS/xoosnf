process XOOS_PARSE_TUMOR_FRACTION_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.1.0'
    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(results_dir)

    output:
    tuple val(meta), path("${meta.id}.tumor_fraction.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-tumor-fraction-metrics \\
        --sample-name '${meta.id}' \\
        --input '${results_dir}/tumor_fraction.results.tsv' \\
        --output '${meta.id}.tumor_fraction.json'
    """

    stub:
    """
    echo '[]' > '${meta.id}.tumor_fraction.json'
    """
}
