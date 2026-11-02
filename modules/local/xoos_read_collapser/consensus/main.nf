process XOOS_READ_COLLAPSER_CONSENSUS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/read_collapser:1.2.0'

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
    def bed_arg = bed ? "--bed ${bed}" : ""
    """
    read_collapser consensus \\
        --bam ${bam} \\
        --output-dir '${meta.id}' \\
        ${bed_arg} \\
        --threads ${task.cpus} \\
        ${args}
    # Collect the read_collapser metric TSVs into a dedicated directory so they
    # can be published without the accompanying consensus FASTQs. The three TSVs
    # are required: a missing file means read_collapser changed its output
    # contract, so fail the task loudly rather than publish a partial set (the
    # silent-absence failure mode this metrics publishing fixes).
    # read_collapser writes metrics under a 'metrics/' subdirectory; preserve
    # that nesting in the emitted dir so markdup/consensus match the
    # rescue-secondary layout and the tool's own output hierarchy.
    mkdir -p read_collapser/metrics
    for f in summary_stats.tsv cluster_size_distribution_summary.tsv cluster_size_distributions.tsv; do
        mv "${meta.id}/metrics/\$f" read_collapser/metrics/
    done
    """

    stub:
    """
    mkdir -p '${meta.id}'
    touch '${meta.id}/output.R1.fastq.gz'
    touch '${meta.id}/output.R2.fastq.gz'
    mkdir -p read_collapser/metrics
    touch 'read_collapser/metrics/summary_stats.tsv'
    touch 'read_collapser/metrics/cluster_size_distribution_summary.tsv'
    touch 'read_collapser/metrics/cluster_size_distributions.tsv'
    """
}
