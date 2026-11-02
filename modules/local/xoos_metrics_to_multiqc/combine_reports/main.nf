process XOOS_COMBINE_MULTIQC_REPORTS {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.2.0'

    label 'process_single'

    tag "${meta.id}"

    input:
    tuple val(meta), path(report_dirs)

    output:
    tuple val(meta), path("combined_multiqc_reports"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('combine_reports.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p combined_multiqc_reports
    combine_reports.py \\
        --input-dirs ${report_dirs} \\
        --output-dir combined_multiqc_reports \\
        ${args} \\
    """

    stub:
    """
    mkdir -p combined_multiqc_reports
    """
}
