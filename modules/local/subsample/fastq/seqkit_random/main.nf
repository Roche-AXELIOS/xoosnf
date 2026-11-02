process SUBSAMPLE_FASTQ_RANDOM {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::seqkit=2.13.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.13.0--he881be0_0' :
        'quay.io/biocontainers/seqkit:2.13.0--he881be0_0' }"

    input:
    tuple val(meta), val(read_count), path(fastqs, stageAs: 'fastq/*')

    output:
    tuple val(meta), path("subsampled/${meta.id}.subsampled.fastq.gz"), emit: fastq
    tuple val("${task.process}"), val('seqkit'), val("2.13.0"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail
    mkdir -p subsampled

    # The -2/--two-pass flag reads the files twice to keep
    # memory usage low on large FASTQs (requires real files, not stdin).
    find -L 'fastq' -name '*.fastq.gz' | sort > fastq_list.txt
    [[ -s fastq_list.txt ]] || { echo "ERROR: no *.fastq.gz files found in fastq/" >&2; exit 1; }

    seqkit sample2 \\
        -2 \\
        -n ${read_count} \\
        -j ${task.cpus} \\
        ${args} \\
        --infile-list fastq_list.txt \\
        -o subsampled/${meta.id}.subsampled.fastq.gz
    """

    stub:
    """
    mkdir -p subsampled
    echo "" | gzip > subsampled/${meta.id}.subsampled.fastq.gz
    """
}
