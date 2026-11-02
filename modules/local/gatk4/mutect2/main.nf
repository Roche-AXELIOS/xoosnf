process GATK4_MUTECT2 {
    tag "${meta.id}"
    label 'process_low'

    conda false
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.3.0.0--py36hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai), path(gzi)
    tuple val(meta4), path(dict)
    path alleles
    path alleles_tbi
    path germline_resource
    path germline_resource_tbi
    path panel_of_normals
    path panel_of_normals_tbi

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi"), emit: tbi
    tuple val(meta), path("*.stats"), emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), emit: f1r2, optional: true
    tuple val(meta), path("*.realigned.bam"), emit: bam, optional: true
    tuple val(meta), path("*.realigned.bai"), emit: bai, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = input.collect { vcf_ -> "--input ${vcf_}" }.join(' \\\n')
    def interval_command = intervals ? "--intervals ${intervals}" : ""
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    def pon_command = panel_of_normals ? "--panel-of-normals ${panel_of_normals}" : ""
    def gr_command = germline_resource ? "--germline-resource ${germline_resource}" : ""
    def a_command = alleles ? "--alleles ${alleles}" : ""
    """
    gatk ${args2} \\
        Mutect2 \\
        ${inputs} \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        ${pon_command} \\
        ${gr_command} \\
        ${a_command} \\
        ${interval_command} \\
        ${bamout_command} \\
        --tmp-dir . \\
        --native-pair-hmm-threads '${task.cpus}' \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    def stub_realigned_bai = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bai" : ""
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    echo "" | gzip > ${prefix}.f1r2.tar.gz
    ${stub_realigned_bam}
    ${stub_realigned_bai}
    """
}
