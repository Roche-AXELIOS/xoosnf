def expandMultipleOptions(optionName, values) {
    return values.collect { value -> "${optionName} '${value}'" }.join(' \\\n')
}

def alignerSampleTypeOptions(sample_type, aligner) {
    def sample_type_name = sample_type.toLowerCase()
    if (sample_type_name != 'ffpe' && sample_type_name != 'cell-line') {
        error("Unknown sample type: ${sample_type}")
    }
    // Remove hyphen from sample_type for model file naming
    sample_type_name = sample_type_name.replaceAll('-', '')

    // Normalize aligner to model file naming: pb_giraffe/vg_giraffe -> giraffe, pb_fq2bam/bwa -> bwa
    def aligner_name = aligner.toLowerCase()
    if (aligner_name in ['pb_giraffe', 'vg_giraffe']) {
        aligner_name = 'giraffe'
    } else if (aligner_name in ['pb_fq2bam', 'bwa']) {
        aligner_name = 'bwa'
    } else {
        log.warn("Unknown aligner '${aligner}' for SVC model selection, defaulting to BWA model")
        aligner_name = 'bwa'
    }

    return "--sample-type ${sample_type_name} --aligner ${aligner_name}"
}

process XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_TUMOR_NORMAL_WGS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.1.0'

    label 'process_high'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bams), path(bam_indexes), path(vcf), path(vcf_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(model)
    tuple val(meta4), path(popaf), path(popaf_index)
    tuple val(meta6), path(target_regions_bed)
    val(aligner)

    output:
    tuple val(meta), path("${prefix}_svc.vcf.gz"), path("${prefix}_svc.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('filter_variants'), eval("filter_variants --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    filter_variants \\
        tumor-normal-wgs \\
        ${expandMultipleOptions('--bam-input', bams)} \\
        --vcf-input '${vcf}' \\
        ${popaf ? "--pop-af-vcf ${popaf}" : ''} \\
        ${model ? "--model '${model}'" : ''} \\
        ${alignerSampleTypeOptions(meta.sample_type ?: "ffpe", aligner)} \\
        --genome '${reference}' \\
        ${target_regions_bed ? "--target-regions '${target_regions_bed}'" : ''} \\
        --output-dir '.' \\
        --vcf-output '${prefix}_svc.vcf.gz' \\
        --tumor-sample-name '${meta.tumor_id ?: meta.id}' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    gzip -n < /dev/null > '${prefix}_svc.vcf.gz'
    touch '${prefix}_svc.vcf.gz.tbi'
    """
}
