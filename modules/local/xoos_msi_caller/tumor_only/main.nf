process XOOS_MSI_CALLER_TUMOR_ONLY {
    conda false
    container 'ghcr.io/roche-axelios/xoos/msi_caller:1.2.0'

    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bam_index)
    tuple val(meta2), path(background_threshold_json)

    output:
    tuple val(meta), path("msi_caller/${prefix}_msi_results.vcf.gz"), emit: vcf
    tuple val(meta), path("msi_caller/${prefix}_msi_results.json"), emit: json
    tuple val(meta), path("msi_caller/metrics/${prefix}_msi_statistics.tsv"), emit: tsv
    tuple val(meta), path("msi_caller/*"), emit: results_dir
    tuple val("${task.process}"), val('msi_caller'), eval("msi_caller --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    // When supplied, Nextflow stages the override and passes its local path.
    // Otherwise msi_caller uses the resource bundled in its container.
    def background_args = background_threshold_json ? "--background-threshold-json '${background_threshold_json}'" : ''
    """
    mkdir -p "msi_caller"

    msi_caller tumor-only \\
        --tumor-bam '${tumor_bam}' \\
        --sample-name '${meta.tumor_sample_name ?: meta.tumor_id ?: meta.id}' \\
        --output-dir 'msi_caller' \\
        --output-prefix '${prefix}' \\
        ${background_args} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // Surface the constructed ext.args into the stub JSON so nf-test can assert
    // the effective tunable flags (e.g. --min-depth) without a container.
    def args = task.ext.args ?: ''
    """
    mkdir -p "msi_caller"
    touch 'msi_caller/${prefix}_msi_results.vcf.gz'
    touch 'msi_caller/${prefix}_msi_results.vcf.gz.tbi'
    echo '{"sample": "${prefix}", "args": "${args}"}' > 'msi_caller/${prefix}_msi_results.json'
    mkdir -p 'msi_caller/metrics'
    touch 'msi_caller/metrics/${prefix}_msi_statistics.tsv'
    """
}
