process XOOS_PARSE_SOMALIER_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.2.0'

    conda false

    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(pairs_tsv), path(samples_tsv), val(sample_ids)

    output:
    tuple val(meta), path("*.somalier*.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), val("1.0"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Pair reports are keyed on the T/N group id while sample reports are keyed
    // on the run name; distinct suffixes keep the two from colliding when both
    // are staged together by the run-level report generator.
    def scope          = pairs_tsv && samples_tsv ? 'somalier' : (pairs_tsv ? 'somalier_pairs' : 'somalier_samples')
    def prefix         = task.ext.prefix ?: "${meta.id}.${scope}"
    def pairs_arg      = pairs_tsv   ? "--pairs-tsv '${pairs_tsv}'"   : ''
    def samples_arg    = samples_tsv ? "--samples-tsv ${samples_tsv}" : ''
    def sample_ids_arg = sample_ids  ? "--sample-ids ${sample_ids.collect { "'${it}'" }.join(' ')}" : ''
    """
    generate-metrics-report parse-somalier-metrics \\
        --key-name '${meta.id}' \\
        ${pairs_arg} \\
        ${samples_arg} \\
        ${sample_ids_arg} \\
        --output '${prefix}.json'
    """

    stub:
    def scope  = pairs_tsv && samples_tsv ? 'somalier' : (pairs_tsv ? 'somalier_pairs' : 'somalier_samples')
    def prefix = task.ext.prefix ?: "${meta.id}.${scope}"
    """
    echo '[]' > '${prefix}.json'
    """
}
