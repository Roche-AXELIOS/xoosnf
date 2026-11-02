process XOOS_BAM_POST_PROCESS_RECRUIT {
    conda false
    container 'ghcr.io/roche-axelios/xoos/bam_post_process:1.2.0'

    label 'process_medium'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index)

    output:
    tuple val(meta), path("${prefix}.recruit.names"), path("${prefix}.recruit.bed"), path("${prefix}.reads.fastq"), emit: recruited
    tuple val(meta), val(args), emit: command_args
    tuple val("${task.process}"), val('process-mitochondrial-bam'), eval("process-mitochondrial-bam --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    process-mitochondrial-bam recruit \\
        --orig-bam '${bam}' \\
        --nr-positions '/resources/numt_catalog.grch38.nr.tsv' \\
        --ucsc-bed '/resources/nuMtSeq_hg38.ucsc.bed' \\
        --out-names '${prefix}.recruit.names' \\
        --out-regions '${prefix}.recruit.bed' \\
        --out-fastq '${prefix}.reads.fastq' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    touch '${prefix}.recruit.names'
    touch '${prefix}.recruit.bed'
    touch '${prefix}.reads.fastq'
    """
}
