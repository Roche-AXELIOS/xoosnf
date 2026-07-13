process XOOS_ALIGNMENT_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.1.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(read_metrics), path(coverage_metrics), path(accuracy_metrics)

    output:
    tuple val(meta), path("*_alignment_metrics_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('xoos_to_multiqc.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Region and pipeline stage label the read/coverage tab so the report can
    // show all regions/stages together. Derived from meta set by the subworkflow.
    def region = meta.metrics_region ?: ''
    def stage = meta.report_prefix ?: ''
    def region_arg = region ? "--region '${region}'" : ''
    def stage_arg = stage ? "--stage '${stage}'" : ''
    // Output dir is unique per sample+stage+region so conversion outputs from
    // different regions do not overwrite each other before they are combined.
    def label = [meta.id, stage, region].findAll { it }.join('_')
    def prefix = task.ext.prefix ?: label
    """
    mkdir -p '${prefix}_alignment_metrics_multiqc_input'
    xoos_to_multiqc.py \\
        --sample-name '${meta.id}' \\
        ${region_arg} \\
        ${stage_arg} \\
        --accuracy-metrics ${accuracy_metrics}/base_level_accuracy_summary.tsv \\
        --hp-errors ${accuracy_metrics}/hp_errors.tsv \\
        --coverage-metrics ${coverage_metrics}/coverage_stats.tsv \\
        --read-metrics ${read_metrics}/read_metrics_summary.tsv \\
        --read-length-summary ${read_metrics}/read_length_summary.tsv \\
        --output-dir '${prefix}_alignment_metrics_multiqc_input' \\
        ${args}
    """

    stub:
    def region = meta.metrics_region ?: ''
    def stage = meta.report_prefix ?: ''
    def label = [meta.id, stage, region].findAll { it }.join('_')
    def prefix = task.ext.prefix ?: label
    """
    mkdir -p '${prefix}_alignment_metrics_multiqc_input'
    touch '${prefix}_alignment_metrics_multiqc_input/xoos_accuracy_metrics_mqc.json'
    touch '${prefix}_alignment_metrics_multiqc_input/xoos_coverage_metrics_mqc.json'
    touch '${prefix}_alignment_metrics_multiqc_input/xoos_read_metrics_read_counts_mqc.json'
    touch '${prefix}_alignment_metrics_multiqc_input/xoos_read_metrics_base_counts_mqc.json'
    """
}
