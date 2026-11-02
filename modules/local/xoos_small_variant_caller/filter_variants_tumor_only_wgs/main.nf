def expandMultipleOptions(optionName, values) {
    return values.collect { value -> "${optionName} '${value}'" }.join(' \\\n')
}

process XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_TUMOR_ONLY_WGS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.2.0'

    label 'process_high'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bams), path(bam_indexes), path(vcf), path(vcf_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(model)
    tuple val(meta4), path(popaf), path(popaf_index)
    tuple val(meta5), path(pon), path(pon_index)
    tuple val(meta6), path(target_regions_bed)

    output:
    tuple val(meta), path("${prefix}_svc.vcf.gz"), path("${prefix}_svc.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}_svc_metrics.tsv"), emit: metrics
    tuple val("${task.process}"), val('filter_variants'), eval("filter_variants --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    // --model is omitted by default so filter_variants uses the tumor-only-wgs
    // model shipped in the container under /resources/. Provide a model file to
    // override. --pon-vcf (panel of normals) is optional; when omitted the PON
    // features default to 0.
    """
    filter_variants \\
        tumor-only-wgs \\
        ${expandMultipleOptions('--bam-input', bams)} \\
        --vcf-input '${vcf}' \\
        ${popaf ? "--pop-af-vcf ${popaf}" : ''} \\
        ${pon ? "--pon-vcf '${pon}'" : ''} \\
        ${model ? "--model '${model}'" : ''} \\
        --genome '${reference}' \\
        ${target_regions_bed ? "--target-regions '${target_regions_bed}'" : ''} \\
        --output-dir '.' \\
        --vcf-output '${prefix}_svc.vcf.gz' \\
        --output-metrics '${prefix}_svc_metrics.tsv' \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    gzip -n < /dev/null > '${prefix}_svc.vcf.gz'
    touch '${prefix}_svc.vcf.gz.tbi'
    printf 'metric_name\\tvalue\\tpercentage\\tdenominator\\nqc_status\\tNA\\tNA\\tNA\\n' > '${prefix}_svc_metrics.tsv'
    """
}
