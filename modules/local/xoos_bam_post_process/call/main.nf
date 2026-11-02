process XOOS_BAM_POST_PROCESS_CALL {
    conda false
    container 'ghcr.io/roche-axelios/xoos/bam_post_process:1.2.0'

    label 'process_medium'

    tag "${meta.id}"

    input:
    tuple val(meta), path(gaf), path(paths_gaf), path(orig_bam), path(orig_bam_index), path(recruit_regions)

    output:
    tuple val(meta), path("${prefix}.chrM.gaf"), emit: gaf
    tuple val(meta), val(args), emit: command_args
    tuple val("${task.process}"), val('process-mitochondrial-bam'), eval("process-mitochondrial-bam --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    process-mitochondrial-bam call \\
        --gaf '${gaf}' \\
        --paths-gaf '${paths_gaf}' \\
        --orig-bam '${orig_bam}' \\
        --recruit-regions '${recruit_regions}' \\
        --nr-positions '/resources/numt_catalog.grch38.nr.tsv' \\
        --control-regions '/resources/control_regions.grch38.bed' \\
        --out-gaf '${prefix}.chrM.gaf' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    touch '${prefix}.chrM.gaf'
    """
}
