process PARABRICKS_FQ2BAM_SINGLE_END_MULTIPLE_FASTQS {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'

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
    def sample_id = meta.sample_id ?: "${meta.id}"
    def run_type  = meta.run_type
    // Parabricks has no --read-group-pm flag, so the full read group (including
    // the PM platform-model tag) is passed as a single quoted @RG string
    // immediately after the input FASTQ, per the pbrun --in-*-fq syntax. This
    // also keeps a single sample-level read group. \\t renders as a literal
    // tab escape that pbrun/bwa expands into @RG field separators.
    def read_group = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:unknown.${run_type}.${sample_id}\\tPL:SBX\\tPU:unknown\\tPM:AXELIOS_1"

    def known_sites_command = known_sites ? (known_sites instanceof List ? known_sites.collect { knownSite -> "--knownSites ${knownSite}" }.join(' \\\n') : "--knownSites ${known_sites}") : ""
    def known_sites_output_cmd = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def intervals_command = intervals ? (intervals instanceof List ? intervals.collect { interval -> "--interval-file ${interval}" }.join(' \\\n') : "--interval-file ${intervals}") : ""

    """
    INDEX=`find -L ${index} -name "*.amb" | sed 's/\\.amb\$//'`

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
        fq2bam \\
        --ref \$INDEX \\
        --in-se-fq 'input.fastq.gz' '${read_group}' \\
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
