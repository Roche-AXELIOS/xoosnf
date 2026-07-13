process PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS {
    conda false
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1'
        : 'nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1'}"

    // For the high mem configuration Giraffe needs at least:
    // GPU memory: 70GB
    // Host memory: 128GB given 96GB as --memory-limit
    label 'process_high'
    label 'process_high_memory'
    label 'process_gpu'

    tag "${meta.id}"

    input:
    tuple val(meta), path(reads, stageAs: 'fastq/*'), path(graph), path(zipcodes), path(minimizers), path(dist)
    path(ref_paths)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample_id = meta.sample_id
    def ref_paths_arg = ref_paths ? "--ref-paths '${ref_paths}'" : ''
    """
    rm -f input.fastq.gz
    mkfifo input.fastq.gz
    find -L 'fastq' -name '*.fastq.gz' -print0 | xargs -0 -n 10 cat > input.fastq.gz &

    pbrun giraffe \\
        --gbz-name '${graph}' \\
        --minimizer-name '${minimizers}' \\
        --dist-name '${dist}' \\
        --zipcodes-name '${zipcodes}' \\
        --in-se-fq input.fastq.gz \\
        --out-bam '${meta.id}.bam' \\
        --sample '${sample_id}' \\
        --read-group '${sample_id}' \\
        ${ref_paths_arg} \\
        ${args}
    """

    stub:
    """
    touch '${meta.id}.bam'
    touch '${meta.id}.bam.bai'
    """
}
