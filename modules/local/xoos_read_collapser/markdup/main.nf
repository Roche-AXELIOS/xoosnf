process XOOS_READ_COLLAPSER_MARK_DUPLICATES {
    conda false
    container 'ghcr.io/roche-axelios/xoos/read_collapser:1.1.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index)

    output:
    tuple val(meta), path("${meta.id}_markdup.bam"), path("${meta.id}_markdup.bam.bai"), path("read_collapser"), emit: markdup_bams
    tuple val(meta), path("read_collapser"), emit: metrics_dir
    tuple val("${task.process}"), val('read_collapser'), eval("read_collapser --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    read_collapser markdup \\
        --bam-input ${bam} \\
        --output-dir rc_work \\
        ${args}
    mv rc_work/output.bam ${meta.id}_markdup.bam
    mv rc_work/output.bam.bai ${meta.id}_markdup.bam.bai
    # Collect only the read_collapser metric TSVs into the published directory so
    # no BAMs or intermediate artifacts leak into the emitted metrics_dir. This
    # mirrors the consensus path's explicit allowlist. The three TSVs are
    # required: a missing file means read_collapser changed its output contract,
    # so fail the task loudly rather than publish a partial set.
    mkdir -p read_collapser
    for f in summary_stats.tsv cluster_size_distribution_summary.tsv cluster_size_distributions.tsv; do
        mv "rc_work/\$f" read_collapser/
    done
    """

    stub:
    """
    touch '${meta.id}_markdup.bam'
    touch '${meta.id}_markdup.bam.bai'
    mkdir -p read_collapser
    touch 'read_collapser/summary_stats.tsv'
    touch 'read_collapser/cluster_size_distribution_summary.tsv'
    touch 'read_collapser/cluster_size_distributions.tsv'
    """
}
