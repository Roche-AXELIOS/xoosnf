process KMC_KMER_COUNT {
    conda false
    container "ghcr.io/roche-axelios/xoosnf/kmc:1.2.0"

    tag "${meta.id}"

    label 'process_high'

    input:
    tuple val(meta), path(fastqs, stageAs: 'fastq/*')
    val read_count

    output:
    tuple val(meta), path("*_personalized_graph.kff"), emit: kff
    tuple val("${task.process}"), val('kmc'), eval("echo \$kmc_version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    find -L ${fastqs} -name '*.fastq.gz' > input_fastq_list.txt
    fastq_downsample \\
        --input-file input_fastq_list.txt \\
        --read-count '${read_count}' \\
        --output-dir 'downsampled_fastq' \\
        --output-file-list 'input_downsample_fastq_list.txt' \\
        --threads ${task.cpus} \\
        ${args2}

    # KMC has a bug when reading from multiple gzipped FASTQ files, so we concatenate them into a single gzipped FASTQ file
    # until this bug is fixed
    cat input_downsample_fastq_list.txt | xargs -n 10 cat > '${meta.id}.downsampled.fastq.gz'

    mkdir -p tmp
    kmc \\
        -t${task.cpus} \\
        ${args} \\
        '${meta.id}.downsampled.fastq.gz' \\
        '${meta.id}_personalized_graph' \\
        tmp

    kmc_version=\$(kmc -v 2>&1 | head -n 1 | grep -oE '[0-9]+([.][0-9]+)+' | head -1)
    """

    stub:
    """
    kmc_version="stub"
    touch '${meta.id}_personalized_graph.kff'
    """
}
