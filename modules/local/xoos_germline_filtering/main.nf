process XOOS_GERMLINE_FILTERING {
    conda false
    container 'ghcr.io/roche-axelios/xoos/germline_filtering:1.2.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf), path(seg)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}_prediction.csv"), emit: prediction
    tuple val(meta), path("metric/${prefix}.tsv"), optional: true, emit: metrics
    tuple val("${task.process}"), val('germline_filtering'), eval("germline_filtering --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    germline_filtering predict \\
        --vcf '${vcf}' \\
        --seg '${seg}' \\
        --output-prefix '${prefix}' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p metric
    gzip -n < /dev/null > '${prefix}.vcf.gz'
    touch '${prefix}.vcf.gz.tbi'
    touch '${prefix}_prediction.csv'
    touch 'metric/${prefix}.tsv'
    """
}
