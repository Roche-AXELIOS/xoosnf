process XOOS_STR_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.2.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(results_files)

    output:
    tuple val(meta), path("*_str_metrics_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('xoos_to_multiqc.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    metrics_file=\$(find -L . -maxdepth 2 -name 'metrics.tsv' | head -1)
    if [ -z "\$metrics_file" ]; then
        echo "ERROR: metrics.tsv not found in staged files" >&2
        exit 1
    fi

    mkdir -p '${prefix}_str_metrics_multiqc_input'
    xoos_to_multiqc.py \\
        --sample-name ${prefix} \\
        --str-metrics "\$metrics_file" \\
        --output-dir '${prefix}_str_metrics_multiqc_input' \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p '${prefix}_str_metrics_multiqc_input'
    touch '${prefix}_str_metrics_multiqc_input/xoos_str_metrics_mqc.json'
    """
}
