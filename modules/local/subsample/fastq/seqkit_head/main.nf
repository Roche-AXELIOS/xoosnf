process SUBSAMPLE_FASTQ {
    tag "${meta.id}"
    label 'process_low'

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
    mkdir -p subsampled

    # Build a sorted file list, then pass it to seqkit via --infile-list.
    # This avoids piping through cat (which causes SIGPIPE when seqkit head
    # exits early) and avoids ARG_MAX limits with many files.
    find -L 'fastq' -name '*.fastq.gz' | sort > fastq_list.txt
    [[ -s fastq_list.txt ]] || { echo "ERROR: no *.fastq.gz files found in fastq/" >&2; exit 1; }
    seqkit head -n ${read_count} -j ${task.cpus} ${args} --infile-list fastq_list.txt -o subsampled/${meta.id}.subsampled.fastq.gz
    """

    stub:
    """
    mkdir -p subsampled
    echo "" | gzip > subsampled/${meta.id}.subsampled.fastq.gz
    """
}
