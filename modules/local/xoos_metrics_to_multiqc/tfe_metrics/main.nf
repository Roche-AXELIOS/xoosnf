process XOOS_TFE_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.1.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(tfe_results_dir)

    output:
    tuple val(meta), path("*_tfe_metrics_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('xoos_to_multiqc.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    mkdir -p '${prefix}_tfe_metrics_multiqc_input'
    xoos_to_multiqc.py \\
        --sample-name ${prefix} \\
        --tfe-results ${tfe_results_dir}/tumor_fraction.results.tsv \\
        --output-dir '${prefix}_tfe_metrics_multiqc_input' \\
        ${args}
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p '${prefix}_tfe_metrics_multiqc_input'
    touch '${prefix}_tfe_metrics_multiqc_input/xoos_tfe_metrics_mqc.json'
    """
}
