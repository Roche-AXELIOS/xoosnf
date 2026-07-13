process VG_HAPLOTYPES {
    conda false
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/vgteam/vg:v1.70.0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    tag "${meta.id}"

    label 'process_high'
    label 'process_high_memory'

    input:
    tuple val(meta), path(kmer_dict)
    tuple val(meta2), path(giraffe_index_dir)

    output:
    tuple val(meta), path("${meta.id}_personalized_graph.gbz"), emit: graphs
    tuple val("${task.process}"), val('vg'), eval("vg version | head -1 | sed 's/vg version //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    REFERENCE_PREFIX=`find -L ${giraffe_index_dir} -name "*.gbz" | head -1 | xargs basename | sed 's/\\.gbz\$//'`

    vg haplotypes \\
        ${args} \\
        -i ${giraffe_index_dir}/\$REFERENCE_PREFIX.hapl \\
        -k ${kmer_dict} \\
        -g ${meta.id}_personalized_graph.gbz \\
        ${giraffe_index_dir}/\$REFERENCE_PREFIX.gbz
    """

    stub:
    """
    touch '${meta.id}_personalized_graph.gbz'
    """
}
