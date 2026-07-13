process SAMBAMBA_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::sambamba=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sambamba:1.0.1--he614052_4' :
        'quay.io/biocontainers/sambamba:1.0.1--he614052_4' }"

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${prefix}.merged.bam"), path("${prefix}.merged.bam.bai"), emit: bam
    tuple val("${task.process}"), val('sambamba'), eval("sambamba --version 2>&1 | head -1 | sed 's/sambamba //'"), topic: versions, emit: versions_sambamba

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sambamba merge \
        ${args} \
        ${prefix}.merged.bam \
        ${bams}

    sambamba index \
        ${args2} \
        ${prefix}.merged.bam \
        ${prefix}.merged.bam.bai
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.bam
    touch ${prefix}.merged.bam.bai
    """
}
