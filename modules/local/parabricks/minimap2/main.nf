process PARABRICKS_MINIMAP2 {
    conda false
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'

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
    def sample_id = meta.sample_id ?: "${meta.id}"
    def run_type  = meta.run_type
    // Parabricks has no --read-group-pm flag, so the full read group (including
    // the PM platform-model tag) is passed as a single quoted @RG string
    // immediately after the input FASTQ, per the pbrun --in-fq syntax. This
    // also keeps a single sample-level read group. \\t renders as a literal
    // tab escape that pbrun expands into @RG field separators.
    def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:unknown.${run_type}.${sample_id}\\tPL:SBX\\tPU:unknown\\tPM:AXELIOS_1"

    def known_sites_command    = known_sites ? (known_sites instanceof List ? known_sites.collect { knownSite -> "--knownSites ${knownSite}" }.join(' \\\n') : "--knownSites ${known_sites}") : ""
    def known_sites_output_cmd = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def intervals_command  = intervals     ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' \\\n') : "--interval-file ${intervals}") : ""

    """
    # Fail fast if no input FASTQs were staged. Without this the writer below
    # would create an empty stream and pbrun would run with no reads.
    if [ -z "\$(find -L 'fastq' -name '*.fastq.gz' -print -quit 2>/dev/null)" ]; then
        echo 'ERROR: no FASTQ files (*.fastq.gz) found under fastq/' >&2
        exit 1
    fi

    # Stream all FASTQ shards into pbrun through a single named pipe (FIFO):
    #   1. Single read group - pbrun emits one @RG per input argument, so a
    #      single concatenated stream yields one sample-level read group
    #      instead of one read group per FASTQ shard.
    #   2. ARG_MAX - passing one FIFO filename avoids the command-line-length
    #      limit that expanding many shard paths onto the command line would hit.
    #   3. Cleanup trap - if pbrun exits before opening the FIFO for reading
    #      (e.g. an argument-parsing error) the background writer would block
    #      forever on open(fifo, O_WRONLY), be orphaned to PID 1, and keep the
    #      container alive until timeout; the EXIT trap kills it so the task
    #      fails fast.
    rm -f input.fastq.gz
    mkfifo input.fastq.gz
    find -L 'fastq' -name '*.fastq.gz' -print0 \\
    | xargs -0 -n 10 cat > input.fastq.gz &
    fifo_writer_pid=\$!
    trap 'kill \$fifo_writer_pid 2>/dev/null || true' EXIT

    pbrun \\
        minimap2 \\
        --ref '${fasta}' \\
        --in-fq input.fastq.gz '${read_group}' \\
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
