process XOOS_SMALL_VARIANT_CALLER_TRAIN_MODEL {
    conda false
    container 'ghcr.io/roche-axelios/xoos/small_variant_caller:1.1.0'

    label 'process_high'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam_features), path(vcf_features), path(truth_vcf), path(truth_vcf_index)

    output:
    tuple val(meta), path("small_variant_caller/${prefix}.model_snv.txt"), emit: model_snv
    tuple val(meta), path("small_variant_caller/${prefix}.model_indel.txt"), emit: model_indel
    tuple val("${task.process}"), val('train_model'), eval("train_model --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "small_variant_caller"

    # Create feature list files
    printf '%s\\n' '${bam_features}' > '${prefix}.bam_features_list.txt'
    printf '%s\\n' '${vcf_features}' > '${prefix}.vcf_features_list.txt'
    printf '%s\\n' '${truth_vcf}' > '${prefix}.vcf_truth_list.txt'

    train_model germline \\
        --positive-bam-features '${prefix}.bam_features_list.txt' \\
        --positive-vcf-features '${prefix}.vcf_features_list.txt' \\
        --truth-vcfs '${prefix}.vcf_truth_list.txt' \\
        --snv-output-file 'small_variant_caller/${prefix}.model_snv.txt' \\
        --indel-output-file 'small_variant_caller/${prefix}.model_indel.txt' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "small_variant_caller"
    touch 'small_variant_caller/${prefix}.model_snv.txt'
    touch 'small_variant_caller/${prefix}.model_indel.txt'
    """
}
