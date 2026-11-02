include { germlineLibraryTypeOption } from '../germline_library_type'

process XOOS_SMALL_VARIANT_CALLER_COMPUTE_VCF_FEATURES {
    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.2.0'

    label 'process_high'

    tag "${meta.id}"

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(target_regions_bed)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(popaf), path(popaf_index)

    output:
    tuple val(meta), path("small_variant_caller/${prefix}.vcf_features.txt"), emit: vcf_features
    tuple val(meta), path("small_variant_caller/${prefix}.vcf.bed"), emit: vcf_bed
    tuple val("${task.process}"), val('compute_vcf_features'), eval("compute_vcf_features --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def target_regions_arg = target_regions_bed ? "--target-regions '${target_regions_bed}'" : ''
    def args = task.ext.args ?: ''
    """
    mkdir -p "small_variant_caller"

    compute_vcf_features germline-wgs \\
        ${germlineLibraryTypeOption(meta.run_type)} \\
        --vcf-input '${vcf}' \\
        --genome '${reference}' \\
        --output-file 'small_variant_caller/${prefix}.vcf_features.txt' \\
        ${target_regions_arg} \\
        --pop-af-vcf '${popaf}' \\
        --output-bed 'small_variant_caller/${prefix}.vcf.bed' \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "small_variant_caller"
    touch 'small_variant_caller/${prefix}.vcf_features.txt'
    touch 'small_variant_caller/${prefix}.vcf.bed'
    """
}
