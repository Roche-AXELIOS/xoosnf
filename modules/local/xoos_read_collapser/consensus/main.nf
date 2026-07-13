process XOOS_READ_COLLAPSER_CONSENSUS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/read_collapser:1.1.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index), path(bed)

    output:
    tuple val(meta), path("${meta.id}/*.fastq.gz"), emit: fastq
    tuple val(meta), path("read_collapser"), emit: metrics_dir
    tuple val("${task.process}"), val('read_collapser'), eval("read_collapser --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bed_arg = bed ? "--bed-input ${bed}" : ""
    """
    read_collapser consensus \\
        --bam-input ${bam} \\
        --output-dir '${meta.id}' \\
        ${bed_arg} \\
        ${args}
    # Collect the read_collapser metric TSVs into a dedicated directory so they
    # can be published without the accompanying consensus FASTQs. The markdup
    # path also emits a clean 'read_collapser' directory; mirror that here so
    # downstream parsing/MultiQC/publishing behave identically across strategies.
    # The three TSVs are required: a missing file means read_collapser changed
    # its output contract, so fail the task loudly rather than publish a partial
    # set (the silent-absence failure mode this metrics publishing fixes).
    mkdir -p read_collapser
    for f in summary_stats.tsv cluster_size_distribution_summary.tsv cluster_size_distributions.tsv; do
        mv "${meta.id}/\$f" read_collapser/
    done
    """

    stub:
    """
    mkdir -p '${meta.id}'
    touch '${meta.id}/output.R1.fastq.gz'
    touch '${meta.id}/output.R2.fastq.gz'
    mkdir -p read_collapser
    touch 'read_collapser/summary_stats.tsv'
    touch 'read_collapser/cluster_size_distribution_summary.tsv'
    touch 'read_collapser/cluster_size_distributions.tsv'
    """
}
