def expandMultipleOptions(optionName, values) {
    return values.collect { value -> "${optionName} '${value}'" }.join(' \\\n')
}

process XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_MITOCHONDRIA {
    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.2.0'

    label 'process_medium'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bams), path(bam_indexes), path(vcf), path(vcf_index)
    tuple val(meta2), path(reference), path(reference_index)

    output:
    tuple val(meta), path("small_variant_caller/${prefix}_mito_svc.vcf.gz"), path("small_variant_caller/${prefix}_mito_svc.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("small_variant_caller/${prefix}_mito_svc_metrics.tsv"), emit: metrics
    tuple val("${task.process}"), val('filter_variants'), eval("filter_variants --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    mkdir -p "small_variant_caller"

    filter_variants \\
        germline-mitochondria \\
        ${expandMultipleOptions('--bam-input', bams)} \\
        --vcf-input '${vcf}' \\
        --genome '${reference}' \\
        --output-dir 'small_variant_caller' \\
        --vcf-output '${prefix}_mito_svc.vcf.gz' \\
        --output-metrics '${prefix}_mito_svc_metrics.tsv' \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p "small_variant_caller"
    echo "" | gzip > 'small_variant_caller/${prefix}_mito_svc.vcf.gz'
    touch 'small_variant_caller/${prefix}_mito_svc.vcf.gz.tbi'
    printf 'metric_name\\tvalue\\tpercentage\\tdenominator\\nqc_status\\tNA\\tNA\\tNA\\n' > 'small_variant_caller/${prefix}_mito_svc_metrics.tsv'
    """
}
