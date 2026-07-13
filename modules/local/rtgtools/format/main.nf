process RTGTOOLS_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dc/dca5ba13b7ec38bf7cacf00a33517b9080067bea638745c05d50a4957c75fc2e/data':
        'community.wave.seqera.io/library/rtg-tools:3.13--3465421f1b0be0ce' }"

    input:
    tuple val(meta), path(input1), path(input2), path(sam_rg)

    output:
    tuple val(meta), path("*.sdf"), emit: sdf
    tuple val("${task.process}"), val('rtgtools'), eval("rtg version 2>/dev/null | sed -n 's/.*Product: RTG Tools //p' | head -n 1"), topic: versions, emit: versions_rtgtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def single = meta.containsKey("single_end") ? meta.single_end : true

    def input = single ? "${input1}" : "--left ${input1} --right ${input2}"
    def rg = sam_rg ? "--sam-rg ${sam_rg}" : ""
    """
    rtg ${args2} format \\
        ${args} \\
        ${rg} \\
        --output ${prefix}.sdf \\
        ${input}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sdf
    """
}
