process PARABRICKS_MINIMAP2 {
    conda false
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1'
        : 'nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1'}"

    label 'process_high'
    label 'process_high_memory'
    label 'process_gpu'

    tag "${meta.id}"

    input:
    tuple val(meta),  path(reads, stageAs: 'fastq/*')
    tuple val(meta2), path(fasta)
    tuple val(meta4), path(intervals)
    tuple val(meta5), path(known_sites)

    output:
    tuple val(meta), path("*.bam"),                   emit: bam,                 optional: true
    tuple val(meta), path("*.bai"),                   emit: bai,                 optional: true
    tuple val(meta), path("*.table"),                 emit: bqsr_table,          optional: true
    tuple val(meta), path("*_qc_metrics"),            emit: qc_metrics,          optional: true
    tuple val(meta), path("*.duplicate-metrics.txt"), emit: duplicate_metrics,   optional: true
    path "compatible_versions.yml",                   emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def known_sites_command    = known_sites ? (known_sites instanceof List ? known_sites.collect { knownSite -> "--knownSites ${knownSite}" }.join(' \\\n') : "--knownSites ${known_sites}") : ""
    def known_sites_output_cmd = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def intervals_command  = intervals     ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' \\\n') : "--interval-file ${intervals}") : ""

    """
    FASTQ_FILES=\$(find -L 'fastq' -name '*.fastq.gz' | tr '\n' ' ')

    pbrun \\
        minimap2 \\
        --ref '${fasta}' \\
        --in-fq \$FASTQ_FILES \\
        --out-bam ${prefix}.bam \\
        ${known_sites_command} \\
        ${known_sites_output_cmd} \\
        ${intervals_command} \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known_sites_output       = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output        = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    ${known_sites_output}
    ${qc_metrics_output}
    ${duplicate_metrics_output}

    # Capture the full version output once and store it in a variable
    pbrun_version_output=\$(pbrun minimap2 --version 2>&1)

    # We handle this different to the other modules because minimap does not begin with an Uppercase letter

    # Generate compatible_versions.yml
    cat <<EOF > compatible_versions.yml
    "${task.process}":
        pbrun_version: \$(echo "\$pbrun_version_output" | grep "pbrun:" | awk '{print \$2}')
        compatible_with:
        \$(echo "\$pbrun_version_output" | tr '\\t' ' ' | awk -F':' '/Compatible With:/,/^---/ { if (\$0 !~ /Compatible With:/ && \$0 !~ /^---\$/ && index(\$0,":")>0) { key=\$1; val=\$2; gsub(/^[ ]+|[ ]+\$/, "", key); gsub(/^[ ]+|[ ]+\$/, "", val); printf "  %s: %s\\n", key, val } }')
    EOF
    """
}
