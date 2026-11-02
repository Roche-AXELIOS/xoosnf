process XOOS_SOMALIER_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.2.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(pairs_tsv), path(samples_tsv)

    output:
    tuple val(meta), path("*_somalier_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), val('1.0'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    def pairs_arg = pairs_tsv ? "--somalier-pairs-tsv ${pairs_tsv}" : ''
    def samples_arg = samples_tsv ? "--somalier-samples-tsv ${samples_tsv}" : ''
    """
    mkdir -p '${prefix}_somalier_multiqc_input'
    xoos_to_multiqc.py \\
        --sample-name ${prefix} \\
        ${pairs_arg} \\
        ${samples_arg} \\
        --output-dir '${prefix}_somalier_multiqc_input' \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p '${prefix}_somalier_multiqc_input'
    touch '${prefix}_somalier_multiqc_input/xoos_somalier_pair_relatedness_mqc.json'
    touch '${prefix}_somalier_multiqc_input/xoos_somalier_sample_stats_mqc.json'
    """
}
