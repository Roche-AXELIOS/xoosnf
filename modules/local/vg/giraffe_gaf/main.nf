process VG_GIRAFFE_GAF {
    conda false
    container 'quay.io/vgteam/vg:v1.70.0'

    label 'process_high'
    label 'process_high_memory'

    tag "${meta.id}"

    input:
    tuple val(meta), path(fastq), path(graph), path(zipcodes), path(minimizers), path(dist)

    output:
    tuple val(meta), path("${prefix}.all.gaf"), emit: gaf
    tuple val(meta), val(args), emit: command_args
    tuple val("${task.process}"), val('vg'), eval("vg version | head -1 | sed 's/vg version //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    vg giraffe \\
        --gbz-name '${graph}' \\
        --minimizer-name '${minimizers}' \\
        --zipcode-name '${zipcodes}' \\
        --dist-name '${dist}' \\
        --fastq-in '${fastq}' \\
        --output-format gaf \\
        ${args} \\
        > '${prefix}.all.gaf'

    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    vg version | head -1 | sed 's/vg version //'
    touch '${prefix}.all.gaf'
    """
}
