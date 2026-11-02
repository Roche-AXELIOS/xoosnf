process VG_SURJECT {
    conda false
    container 'quay.io/vgteam/vg:v1.70.0'

    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(gaf), path(graph)
    val(chromosome)

    output:
    tuple val(meta), path("*.unsorted.bam"), emit: bam
    tuple val("${task.process}"), val('vg'), eval("echo \$vg_version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    vg surject \\
        -x '${graph}' \\
        ${args} \\
        '${gaf}' \\
        > '${prefix}.${chromosome}.unsorted.bam'

    vg_version=\$(vg version | head -1 | sed 's/vg version //')
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    vg_version="stub"
    touch '${prefix}.${chromosome}.unsorted.bam'
    """
}
