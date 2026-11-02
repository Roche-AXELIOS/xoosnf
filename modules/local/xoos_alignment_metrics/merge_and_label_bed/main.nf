process XOOS_MERGE_AND_LABEL_BED {
    conda false
    container 'ghcr.io/roche-axelios/xoos/alignment_metrics:1.2.0'

    label 'process_single'
    tag "${meta.id}"

    input:
    // `labels` and `beds` are positionally paired lists: beds[i] is assigned labels[i].
    // Each BED may be plain or gzip/bgzip compressed (named with a `.gz` extension).
    tuple val(meta), val(labels), path(beds, stageAs: '?/*')

    output:
    tuple val(meta), path("*.labeled.bed"), emit: labeled_bed
    tuple val("${task.process}"), val('alignment_metrics'), eval('merge_and_label_bed --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    def safeLabelPattern = /^[A-Za-z0-9._-]+$/
    def invalidLabels = labels.findAll { label -> !(label.toString() ==~ safeLabelPattern) }
    if (invalidLabels) {
        error("Invalid region label(s): ${invalidLabels.join(', ')}. Labels must match ${safeLabelPattern}.")
    }
    // If only a single BED is provided, wrap it in a list to avoid Nextflow
    // incorrectly interpreted as an iterable of path segments.
    def bedPaths = beds instanceof Collection ? beds : [beds]
    if (bedPaths.size() != labels.size()) {
        error("BED/label count mismatch: received ${bedPaths.size()} BED file(s) and ${labels.size()} label(s).")
    }
    // Pair each BED with its label by position, once per --bed/--label.
    def bed_label_args = [bedPaths, labels]
        .transpose()
        .collect { bed, label -> "--bed '${bed}' --label '${label}'" }
        .join(' \\\n        ')
    """
    merge_and_label_bed \\
        ${bed_label_args} \\
        --output '${prefix}.labeled.bed' \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    install -D /dev/null '${prefix}.labeled.bed'
    """
}
