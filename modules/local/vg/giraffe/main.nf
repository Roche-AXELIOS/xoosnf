process VG_GIRAFFE {
    conda false
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/vgteam/vg:v1.70.0'
        : 'quay.io/vgteam/vg:v1.70.0'}"

    label 'process_high'
    label 'process_high_memory'

    tag "${meta.id}"

    input:
    tuple val(meta), path(fastqs, stageAs: 'fastq/*'), path(graph), path(zipcodes), path(minimizers), path(dist)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val("${task.process}"), val('vg'), eval("vg version | head -1 | sed 's/vg version //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def sample_id = meta.sample_id
    """
    mkdir tmp
    find -L 'fastq' -name '*.fastq.gz' -print0 \\
        | xargs -0 -n 10 zcat \\
        | vg giraffe \\
        --gbz-name '${graph}' \\
        --minimizer-name '${minimizers}' \\
        --zipcode-name '${zipcodes}' \\
        --dist-name '${dist}' \\
        --fastq-in - \\
        --sample '${sample_id}' \\
        ${args} \\
        | samtools sort \\
        ${args2} \\
        -o '${meta.id}.bam##idx##${meta.id}.bam.bai' \\
        -T 'tmp/${meta.id}_sort_tmp'
    """

    stub:
    """
    touch '${meta.id}.bam'
    touch '${meta.id}.bam.bai'
    """
}
