process GATK4_MERGEMUTECTSTATS {
    tag "${meta.id}"
    label 'process_low'

    conda false
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(stats, stageAs: "input_stats/*")

    output:
    tuple val(meta), path("*.stats"), emit: stats
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = stats.collect { s -> "--stats ${s}" }.join(' \\\n')
    """
    gatk ${args} \\
        MergeMutectStats \\
        ${input_list} \\
        --output ${prefix}.vcf.gz.stats \\
        --tmp-dir .
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz.stats
    """
}
