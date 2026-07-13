process PARABRICKS_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    tuple val(ref_meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf, optional: true
    tuple val(meta), path("*.g.vcf.gz"), emit: gvcf, optional: true
    tuple val(meta), path("*_hc*.bam"), path("*_hc*.bam.bai"), emit: bam
    path "compatible_versions.yml", emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit(1, "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def output_file = args.contains("--gvcf") ? "${prefix}.g.vcf.gz" : "${prefix}.vcf"
    def intervals_command = intervals ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' \\\n') : "--interval-file ${intervals}") : ""

    """
    pbrun \\
        haplotypecaller \\
        --ref ${fasta} \\
        --in-bam ${input} \\
        --out-variants ${output_file} \\
        --htvc-bam-output ${prefix}_hc.bam \\
        ${intervals_command} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_cmd = args.contains("--gvcf") ? "echo '' | gzip > ${prefix}.g.vcf.gz" : "touch ${prefix}.vcf"
    def gpu_count = task.accelerator ? task.accelerator.request : 1
    """
    ${output_cmd}
    if [ "$gpu_count" -eq 1 ]; then
        touch ${prefix}_hc.bam
        touch ${prefix}_hc.bam.bai
    else
        for i in \$(seq 0 \$(($gpu_count - 1))); do
            touch ${prefix}_hc_\${i}.bam
            touch ${prefix}_hc_\${i}.bam.bai
        done
    fi

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun haplotypecaller --version 2>&1)

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | awk '/Compatible With:/,/^---/{ if (\$1 ~ /^[A-Z]/ && \$1 != "Compatible" && \$1 != "---") { printf "  %s: %s\\n", \$1, \$2 } }')
    EOF
    """
}
