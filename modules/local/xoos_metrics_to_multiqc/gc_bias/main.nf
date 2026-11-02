process XOOS_GC_BIAS_METRICS_TO_MULTIQC {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/xoos-to-multiqc:1.2.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(gc_bias)

    output:
    tuple val(meta), path("*_gc_bias_multiqc_input"), emit: multiqc_report_dir
    tuple val("${task.process}"), val('xoos_to_multiqc'), eval('xoos_to_multiqc.py --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def stage = meta.report_prefix ?: ''
    def stage_arg = stage ? "--stage '${stage}'" : ''
    // Keep stage-specific outputs distinct when primary and downsampled
    // GC-bias metrics for a sample are combined in one MultiQC report.
    def label = [meta.id, stage].findAll { it }.join('_')
    def prefix = task.ext.prefix ?: label
    """
    mkdir -p '${prefix}_gc_bias_multiqc_input'
    xoos_to_multiqc.py \\
        --sample-name '${meta.id}' \\
        ${stage_arg} \\
        --gc-bias '${gc_bias}/gc_bias.tsv' \\
        --gc-bias-summary '${gc_bias}/gc_bias_summary.tsv' \\
        --output-dir '${prefix}_gc_bias_multiqc_input' \\
        ${args}
    """

    stub:
    def stage = meta.report_prefix ?: ''
    def label = [meta.id, stage].findAll { it }.join('_')
    def prefix = task.ext.prefix ?: label
    """
    install -D /dev/null '${prefix}_gc_bias_multiqc_input/xoos_gc_bias_mqc.json'
    """
}
