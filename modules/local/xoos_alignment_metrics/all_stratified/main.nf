process XOOS_ALIGNMENT_METRICS_ALL_STRATIFIED {
    conda false
    container 'ghcr.io/roche-axelios/xoos/alignment_metrics:1.2.0'

    label 'process_high'
    tag "${meta1.id}"

    input:
    tuple val(meta1), path(bam), path(bam_index)
    tuple val(meta2), path(fasta), path(fai)
    tuple val(meta3), path(bed)

    output:
    // With `--stratify-by-label`, metrics are written to per-label subdirectories:
    //   metrics/<label>/{accuracy,coverage,read_metrics}/*.tsv
    // The whole tree is emitted as a single directory; the calling subworkflow
    // fans it out into per-label metric tuples using the emitted subdirectory names.
    tuple val(meta1), path("metrics"), emit: metrics_dir
    tuple val("${task.process}"), val('alignment_metrics'), eval('alignment_metrics --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reference_opt = fasta ? "--reference ${fasta}" : ""
    def bed_opt = bed ? "--bed ${bed}" : ""
    def args = task.ext.args ?: ''
    """
    mkdir -p metrics
    alignment_metrics all \\
        --bam '${bam}' \\
        ${bed_opt} \\
        ${reference_opt} \\
        --output-dir metrics \\
        --stratify-by-label \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    def labels = meta3.region_labels ?: ['stub_label']
    def make_dirs = labels.collect { label ->
        [
            "install -D /dev/null 'metrics/${label}/coverage/coverage_stats.tsv'",
            "install -D /dev/null 'metrics/${label}/coverage/coverage_distribution_summary.tsv'",
            "install -D /dev/null 'metrics/${label}/accuracy/base_level_accuracy_summary.tsv'",
            "install -D /dev/null 'metrics/${label}/accuracy/hp_errors.tsv'",
            "install -D /dev/null 'metrics/${label}/read_metrics/read_metrics_summary.tsv'",
            "install -D /dev/null 'metrics/${label}/read_metrics/read_length_summary.tsv'",
        ].join('\n    ')
    }.join('\n    ')
    """
    ${make_dirs}
    """
}
