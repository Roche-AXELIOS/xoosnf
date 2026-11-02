process PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS {
    conda false
    container 'nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1'

    // For the high mem configuration Giraffe needs at least:
    // GPU memory: 70GB
    // Host memory: 128GB given 96GB as --memory-limit
    label 'process_high'
    label 'process_high_memory'
    label 'process_gpu'

    tag "${meta.id}"

    input:
    tuple val(meta), path(reads, stageAs: 'fastq/*'), path(graph), path(zipcodes), path(minimizers), path(dist)
    path(ref_paths)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val(meta), val(args), emit: command_args
    tuple val("${task.process}"), val('parabricks'), eval("pbrun version | grep -m1 '^pbrun:' | sed 's/^pbrun:[[:space:]]*//'"), topic: versions, emit: versions_parabricks

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    def sample_id = meta.sample_id ?: "${meta.id}"
    def ref_paths_arg = ref_paths ? "--ref-paths '${ref_paths}'" : ''
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
    find -L 'fastq' -name '*.fastq.gz' -print0 | xargs -0 -n 10 cat > input.fastq.gz &
    fifo_writer_pid=\$!
    trap 'kill \$fifo_writer_pid 2>/dev/null || true' EXIT

    pbrun giraffe \\
        --gbz-name '${graph}' \\
        --minimizer-name '${minimizers}' \\
        --dist-name '${dist}' \\
        --zipcodes-name '${zipcodes}' \\
        --in-se-fq input.fastq.gz \\
        --out-bam '${meta.id}.bam' \\
        --sample '${sample_id}	PM:AXELIOS_1' \\
        --read-group '${sample_id}' \\
        ${ref_paths_arg} \\
        ${args}
    """

    stub:
    args = task.ext.args ?: ''
    """
    touch '${meta.id}.bam'
    touch '${meta.id}.bam.bai'
    """
}
