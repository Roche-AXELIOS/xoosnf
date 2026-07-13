process PARABRICKS_MUTECTCALLER {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(intervals)
    tuple val(ref_meta), path(fasta)
    path panel_of_normals
    path panel_of_normals_index

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"),       emit: vcf
    tuple val(meta), path("*.vcf.gz.stats"), emit: stats
    tuple val(meta), path("*_mutect.bam"), path("*_mutect.bam.bai"), emit: bam
    path "compatible_versions.yml",          emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def intervals_command  = intervals     ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' \\\n') : "--interval-file ${intervals}") : ""
    def prepon_command = panel_of_normals ? "cp -L ${panel_of_normals_index} `readlink -f ${panel_of_normals}`.tbi && pbrun prepon --in-pon-file ${panel_of_normals}" : ""
    def postpon_command = panel_of_normals ? "pbrun postpon --in-vcf ${prefix}.vcf.gz --in-pon-file ${panel_of_normals} --out-vcf ${prefix}_annotated.vcf.gz" : ""

    """
    # if panel of normals specified, run prepon
    ${prepon_command}

    pbrun \\
        mutectcaller \\
        --ref ${fasta} \\
        --in-tumor-bam ${tumor_bam} \\
        --tumor-name ${meta.tumor_id} \\
        --out-vcf ${prefix}.vcf.gz \\
        --mutect-bam-output ${prefix}_mutect.bam \\
        ${intervals_command} \\
        ${args}

    # if panel of normals specified, run postpon
    ${postpon_command}

    # reindex the final VCF, because the one produced by mutectcaller is corrupt
    bcftools index ${args2} ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def postpon_command = panel_of_normals ? "echo '' | gzip > ${prefix}_annotated.vcf.gz" : ""
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    touch ${prefix}_mutect.bam
    touch ${prefix}_mutect.bam.bai
    ${postpon_command}

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun mutectcaller --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF
    """
}
