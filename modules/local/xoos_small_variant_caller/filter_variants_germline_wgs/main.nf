def expandMultipleOptions(optionName, values) {
    return values.collect { value -> "${optionName} '${value}'" }.join(' \\\n')
}

def selectWorkflow(use_multisample_model) {
    if (use_multisample_model) {
        return "germline-multi-sample"
    } else {
        return "germline"
    }
}

def alignerRunTypeOptions(run_type, aligner) {
    def run_type_name = run_type.toLowerCase()
    if (run_type_name != 'sbx-d' && run_type_name != 'sbx-fast') {
        error("Unknown run type: ${run_type}")
    }
    // Remove hyphen from runType for model file naming
    run_type_name = run_type_name.replaceAll('-', '')

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

    return "--run-type ${run_type_name} --aligner ${aligner_name}"
}

process XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.1.0'

    label 'process_high'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bams), path(bam_indexes), path(vcf), path(vcf_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(model_snv), path(model_indel)
    tuple val(meta4), path(popaf), path(popaf_index)
    tuple val(meta5), path(par_bed_x), path(par_bed_y)
    tuple val(meta6), path(target_regions_bed)
    val(use_multisample_model)
    val(aligner)

    output:
    tuple val(meta), path("small_variant_caller/${prefix}_svc.vcf.gz"), path("small_variant_caller/${prefix}_svc.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('filter_variants'), eval("filter_variants --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    """
    mkdir -p "small_variant_caller"

    filter_variants \\
        ${selectWorkflow(use_multisample_model)} \\
        ${expandMultipleOptions('--bam-input', bams)} \\
        --vcf-input '${vcf}' \\
        ${popaf ? "--pop-af-vcf ${popaf}" : ''} \\
        ${model_snv ? "--snv-model '${model_snv}'" : ''} \\
        ${model_indel ? "--indel-model '${model_indel}'" : ''} \\
        ${alignerRunTypeOptions(meta.run_type, aligner)} \\
        --genome '${reference}' \\
        ${par_bed_x ? "--par-bed-x '${par_bed_x}'" : ''} \\
        ${par_bed_y ? "--par-bed-y '${par_bed_y}'" : ''} \\
        ${target_regions_bed ? "--target-regions '${target_regions_bed}'" : ''} \\
        --output-dir 'small_variant_caller' \\
        --vcf-output '${prefix}_svc.vcf.gz' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p "small_variant_caller"
    touch 'small_variant_caller/${prefix}_svc.vcf.gz'
    touch 'small_variant_caller/${prefix}_svc.vcf.gz.tbi'
    """
}
