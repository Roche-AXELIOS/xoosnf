process PARABRICKS_FQ2BAM_SINGLE_END_MULTIPLE_FASTQS {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1'
        : 'nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1'}"

    input:
    tuple val(meta), path(reads, stageAs: "fastq/*")
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    tuple val(meta4), path(intervals)
    tuple val(meta5), path(known_sites)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.table"), emit: bqsr_table, optional: true
    tuple val(meta), path("*_qc_metrics"), emit: qc_metrics, optional: true
    tuple val(meta), path("*.duplicate-metrics.txt"), emit: duplicate_metrics, optional: true
    path "compatible_versions.yml", emit: compatible_versions, optional: true
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_id = meta.sample_id

    def known_sites_command = known_sites ? (known_sites instanceof List ? known_sites.collect { knownSite -> "--knownSites ${knownSite}" }.join(' \\\n') : "--knownSites ${known_sites}") : ""
    def known_sites_output_cmd = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def intervals_command = intervals ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' \\\n') : "--interval-file ${intervals}") : ""

    """
    INDEX=`find -L ${index} -name "*.amb" | sed 's/\\.amb\$//'`
    rm -f input.fastq.gz
    mkfifo input.fastq.gz
    find -L 'fastq' -name '*.fastq.gz' -print0 \\
    | xargs -0 -n 10 cat > input.fastq.gz &

    pbrun \\
        fq2bam \\
        --ref \$INDEX \\
        --in-se-fq 'input.fastq.gz' \\
        --read-group-sm ${sample_id} \\
        --read-group-id-prefix ${sample_id} \\
        --out-bam ${prefix}.bam \\
        ${known_sites_command} \\
        ${known_sites_output_cmd} \\
        ${intervals_command} \\
        ${args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known_sites_output = known_sites ? "touch ${prefix}.table" : ""
    def qc_metrics_output = args.contains("--out-qc-metrics-dir") ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = args.contains("--out-duplicate-metrics") ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    ${known_sites_output}
    ${qc_metrics_output}
    ${duplicate_metrics_output}

    # Capture once and build single-line compatible_with (spaces only, no tabs)
    pbrun_version_output=\$(pbrun fq2bam --version 2>&1)

    # Because of a space between BWA and mem in the version output this is handled different to the other modules
    compat_line=\$(echo "\$pbrun_version_output" | awk -F':' '
        /Compatible With:/ {on=1; next}
        /^---/ {on=0}
        on && /:/ {
            key=\$1; val=\$2
            gsub(/[ \\t]+/, " ", key); gsub(/^[ \\t]+|[ \\t]+\$/, "", key)
            gsub(/[ \\t]+/, " ", val); gsub(/^[ \\t]+|[ \\t]+\$/, "", val)
            a[++i]=key ": " val
        }
        END { for (j=1;j<=i;j++) printf "%s%s", (j>1?", ":""), a[j] }
    ')

    cat <<EOF > compatible_versions.yml
    "${task.process}":
    pbrun_version: \$(echo "\$pbrun_version_output" | awk '/^pbrun:/ {print \$2; exit}')
    compatible_with: "\$compat_line"
    EOF
    """
}
