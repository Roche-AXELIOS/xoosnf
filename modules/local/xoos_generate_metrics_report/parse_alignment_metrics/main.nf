process XOOS_PARSE_ALIGNMENT_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'
    conda false

    label 'process_single'
    tag "${meta.id}, ${source}"

    input:
    tuple val(meta), path(read_metrics_dir), path(coverage_metrics_dir), path(accuracy_metrics_dir)
    val(source)
    val(report_label)

    output:
    tuple val(meta), path("${meta.id}.alignment_metrics_${report_label}.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-alignment-metrics \\
        --sample-name '${meta.id}' \\
        --source 'alignment_metrics/${source}' \\
        --input '${coverage_metrics_dir}/coverage_stats.tsv' \\
        --input '${coverage_metrics_dir}/coverage_distribution_summary.tsv' \\
        --input '${accuracy_metrics_dir}/base_level_accuracy_summary.tsv' \\
        --input '${read_metrics_dir}/read_metrics_summary.tsv' \\
        --input '${read_metrics_dir}/read_length_summary.tsv' \\
        --input '${read_metrics_dir}/read_length_excluding_soft_clips_summary.tsv' \\
        --input '${coverage_metrics_dir}/coverage_uniformity_summary.tsv' \\
        --input '${read_metrics_dir}/read_metrics_te.tsv' \\
        --input-hp-errors '${accuracy_metrics_dir}/hp_errors.tsv' \\
        --output '${meta.id}.alignment_metrics_${report_label}.json'
    """

    stub:
    """
    echo '{}' > '${meta.id}.alignment_metrics_${report_label}.json'
    """
}
