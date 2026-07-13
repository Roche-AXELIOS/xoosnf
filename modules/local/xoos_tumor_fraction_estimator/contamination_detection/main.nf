process XOOS_TUMOR_FRACTION_ESTIMATOR_CONTAMINATION_DETECTION {
    tag "${meta.id}"
    label 'process_high'

    conda false
    container 'ghcr.io/roche-axelios/xoos/tumor_fraction_estimator:1.1.0'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_ref), path(reference), path(fai)
    path(vcf)
    path(included_regions)

    output:
    tuple val(meta), path("contamination_detection"), emit: results
    tuple val("${task.process}"), val('tumor_fraction_estimator'), eval('tumor_fraction_estimator --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def vcf_opt      = vcf ? "--vcf '${vcf}'" : ""
    def included_opt = included_regions ? "--included-regions-bed '${included_regions}'" : "--included-regions-bed '/resources/GRCh38_notinalldifficultregions.noY.bed'"

    """
    mkdir -p 'contamination_detection'

    tumor_fraction_estimator contamination \\
        --bam '${bam}' \\
        ${vcf_opt} \\
        --reference '${reference}' \\
        --output-dir 'contamination_detection' \\
        ${included_opt} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    mkdir -p 'contamination_detection'
    touch 'contamination_detection/contamination.results.tsv'
    touch 'contamination_detection/contamination.counts.tsv'
    """
}
