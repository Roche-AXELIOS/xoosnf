process VG_GIRAFFE {
    conda false
    container 'quay.io/vgteam/vg:v1.70.0'

    label 'process_high'
    label 'process_high_memory'

    tag "${meta.id}"

    input:
    tuple val(meta), path(fastqs, stageAs: 'fastq/*'), path(graph), path(zipcodes), path(minimizers), path(dist)
    path(ref_paths)

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    tuple val(meta), val(args), emit: command_args
    tuple val("${task.process}"), val('vg'), eval("echo \$vg_version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def ref_paths_arg = ref_paths ? "--ref-paths '${ref_paths}'" : ''
    """
    mkdir tmp

    # Fail fast if no input FASTQs were staged. Without this, xargs would run
    # zcat with no file arguments, which reads from stdin and hangs the task.
    if [ -z "\$(find -L 'fastq' -name '*.fastq.gz' -print -quit 2>/dev/null)" ]; then
        echo 'ERROR: no FASTQ files (*.fastq.gz) found under fastq/' >&2
        exit 1
    fi

    find -L 'fastq' -name '*.fastq.gz' -print0 \\
        | xargs -0 -n 10 zcat \\
        | vg giraffe \\
        --gbz-name '${graph}' \\
        --minimizer-name '${minimizers}' \\
        --zipcode-name '${zipcodes}' \\
        --dist-name '${dist}' \\
        --fastq-in - \\
        ${ref_paths_arg} \\
        -t ${task.cpus} \\
        ${args} \\
        | samtools sort \\
        -@ ${task.cpus} \\
        ${args2} \\
        -o '${meta.id}.bam##idx##${meta.id}.bam.bai' \\
        -T 'tmp/${meta.id}_sort_tmp'

    vg_version=\$(vg version | head -1 | sed 's/vg version //')
    """

    stub:
    args = task.ext.args ?: ''
    """
    vg_version="stub"
    touch '${meta.id}.bam'
    touch '${meta.id}.bam.bai'
    """
}
