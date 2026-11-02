process SAMTOOLS_BAM_SAMPLE_NAME {
    tag "${meta.id}:${role}"
    label 'process_single'

    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'

    input:
    tuple val(sample_key), val(meta), val(role), path(bam, stageAs: 'input.bam')

    output:
    tuple val(sample_key), val(meta), val(role), path('sample_name.txt'), emit: sample_name
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    samtools view -H input.bam \\
        | awk -F '\\t' '\$1 == "@RG" { for (i = 2; i <= NF; i++) if (\$i ~ /^SM:/) { sub(/^SM:/, "", \$i); print \$i; exit } }' \\
        > sample_name.txt

    if [ ! -s sample_name.txt ]; then
        echo "ERROR: no SM tag found in the @RG header." >&2
        exit 1
    fi

    # The resolved SM value is propagated into downstream command arguments.
    # Restrict it to shell-safe sample-name characters before emitting it so
    # whitespace and shell metacharacters cannot change argument parsing.
    if ! grep -Eq '^[A-Za-z0-9._-]+\$' sample_name.txt; then
        echo "ERROR: unsupported characters in the SM tag; expected only letters, numbers, '.', '_' or '-'." >&2
        exit 1
    fi
    """

    stub:
    def stub_sample_name = role == 'tumor' ? (meta.tumor_id ?: meta.id) : (meta.normal_id ?: meta.id)
    """
    printf '%s\\n' '${stub_sample_name}' > sample_name.txt
    """
}
