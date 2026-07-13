process FAI_TO_BED {
    tag "$meta.id"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0--2' }"

    input:
    tuple val(meta), path(fasta_fai)

    output:
    tuple val(meta), path("*.bed") , emit: bed
    tuple val("${task.process}"), val('awk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_awk

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'BEGIN {SEP="\t"}; {print \$1 SEP "0" SEP \$2}' ${fasta_fai} > ${fasta_fai.baseName}.bed
    """

    stub:
    """
    touch chr.bed
    """
}
