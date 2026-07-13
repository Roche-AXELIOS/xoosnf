process XOOS_TUMOR_FRACTION_ESTIMATOR_TUMOR_FRACTION {
    tag "${meta.id}"
    label 'process_high'

    conda false
    container 'ghcr.io/roche-axelios/xoos/tumor_fraction_estimator:1.1.0'

    input:
    tuple val(meta), path(bam), path(bai), path(vcf)
    tuple val(meta_ref), path(reference), path(fai)
    tuple val(meta_target), path(included_regions)

    output:
    tuple val(meta), path("tumor_fraction_estimation"), emit: results
    tuple val("${task.process}"), val('tumor_fraction_estimator'), eval('tumor_fraction_estimator --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def included_opt = included_regions ? "--included-regions-bed '${included_regions}'" : "--included-regions-bed '/resources/GRCh38_notinalldifficultregions.noY.bed'"

    """
    mkdir -p 'tumor_fraction_estimation'

    tumor_fraction_estimator tumor-fraction \\
        --bam '${bam}' \\
        --vcf '${vcf}' \\
        --reference '${reference}' \\
        --output-dir 'tumor_fraction_estimation' \\
        ${included_opt} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    mkdir -p tumor_fraction_estimation
    touch "tumor_fraction_estimation/tumor_fraction.results.tsv"
    touch "tumor_fraction_estimation/tumor_fraction.counts.tsv"
    """
}
