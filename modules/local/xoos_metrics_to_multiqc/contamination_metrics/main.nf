process XOOS_CONTAMINATION_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.2.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(contamination_results_dir)

    output:
    tuple val(meta), path("*_contamination_metrics_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('xoos_to_multiqc.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    mkdir -p '${prefix}_contamination_metrics_multiqc_input'
    xoos_to_multiqc.py \\
        --sample-name ${prefix} \\
        --tfe-contamination-results ${contamination_results_dir}/contamination_estimation_results.tsv \\
        --output-dir '${prefix}_contamination_metrics_multiqc_input' \\
        ${args}
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p '${prefix}_contamination_metrics_multiqc_input'
    touch '${prefix}_contamination_metrics_multiqc_input/xoos_tfe_contamination_metrics_mqc.json'
    """
}
