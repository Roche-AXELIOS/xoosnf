process XOOS_DEMUX_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.2.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(run_stats), path(sample_stats)

    output:
    tuple val(meta), path("*_demux_metrics_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('xoos_to_multiqc.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    mkdir -p '${prefix}_demux_metrics_multiqc_input'
    xoos_to_multiqc.py \\
        --demux-run-stats ${run_stats} \\
        --demux-sample-stats ${sample_stats} \\
        --output-dir '${prefix}_demux_metrics_multiqc_input' \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p '${prefix}_demux_metrics_multiqc_input'
    touch '${prefix}_demux_metrics_multiqc_input/xoos_demux_run_metrics_mqc.json'
    touch '${prefix}_demux_metrics_multiqc_input/xoos_demux_sample_counts_mqc.json'
    touch '${prefix}_demux_metrics_multiqc_input/xoos_demux_sample_metrics_mqc.json'
    """
}
