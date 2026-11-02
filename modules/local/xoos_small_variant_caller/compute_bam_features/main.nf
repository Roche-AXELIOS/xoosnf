include { germlineLibraryTypeOption } from '../germline_library_type'

def expandMultipleOptions(optionName, values) {
    return values.collect { value -> "${optionName} '${value}'" }.join(' \\\n')
}

process XOOS_SMALL_VARIANT_CALLER_COMPUTE_BAM_FEATURES {
    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.2.0'

    label 'process_high'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bams), path(bam_indexes), path(target_regions_bed)
    tuple val(meta2), path(reference), path(reference_index)

    output:
    tuple val(meta), path("small_variant_caller/${prefix}.bam_features.txt"), emit: bam_features
    tuple val("${task.process}"), val('compute_bam_features'), eval("compute_bam_features --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def target_regions_arg = target_regions_bed ? "--target-regions '${target_regions_bed}'" : ''
    def args = task.ext.args ?: ''
    """
    mkdir -p "small_variant_caller"

    compute_bam_features germline-wgs \\
        ${germlineLibraryTypeOption(meta.run_type)} \\
        ${expandMultipleOptions('--bam-input', bams)} \\
        --genome '${reference}' \\
        --output-file 'small_variant_caller/${prefix}.bam_features.txt' \\
        ${target_regions_arg} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "small_variant_caller"
    touch 'small_variant_caller/${prefix}.bam_features.txt'
    """
}
