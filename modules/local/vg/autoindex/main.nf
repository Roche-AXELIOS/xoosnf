process VG_AUTOINDEX {
    conda false
    container 'quay.io/vgteam/vg:v1.70.0'

    tag "${meta.id}"

    label 'process_high'
    label 'process_high_memory'

    input:
    tuple val(meta), path("${meta.id}_personalized_graph.gbz")

    output:
    tuple val(meta), path("${meta.id}_personalized_graph.gbz"), path("${meta.id}_personalized_graph.shortread.zipcodes"), path("${meta.id}_personalized_graph.shortread.withzip.min"), path("${meta.id}_personalized_graph.dist"), emit: graphs
    tuple val("${task.process}"), val('vg'), eval("vg version | head -1 | sed 's/vg version //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vg autoindex \\
        -t ${task.cpus} \\
        ${args} \\
        -p '${meta.id}_personalized_graph' \\
        -G '${meta.id}_personalized_graph.gbz'
    """

    stub:
    """
    touch '${meta.id}_personalized_graph.gbz'
    touch '${meta.id}_personalized_graph.shortread.zipcodes'
    touch '${meta.id}_personalized_graph.shortread.withzip.min'
    touch '${meta.id}_personalized_graph.dist'
    """
}
