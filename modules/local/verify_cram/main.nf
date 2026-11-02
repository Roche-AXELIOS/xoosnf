process VERIFY_CRAM {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'

    input:
    tuple val(meta), path(cram), path(crai)

    output:
    tuple val(meta), path(cram), path(crai), emit: cram
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -uo pipefail

    samtools view -H ${cram} > header.sam

    if ! awk -F'\\t' '\$1=="@SQ" && \$0!~/\\tM5:/ { exit 1 }' header.sam; then
        echo "ERROR: ${cram} is missing M5 (MD5 checksum) tags on its @SQ header records." >&2
        exit 1
    fi
    """

    stub:
    // The M5 check is a no-op in stub mode: the emitted cram/crai are the staged
    // inputs, so nothing needs to be created.
    """
    :
    """
}
