process XOOS_TUMOR_FRACTION_ESTIMATOR_TUMOR_FRACTION {
    tag "${meta.id}"
    label 'process_high'

    conda false
    container 'ghcr.io/roche-axelios/xoos/tumor_fraction_estimator:1.2.0'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta_reporter_vcf), path(vcf), path(tbi)
    tuple val(meta_ref), path(reference), path(fai)
    tuple val(meta_target), path(included_regions)
    tuple val(meta_pon), path(panel_of_normals_tsv), path(panel_of_normals_tbi)

    output:
    tuple val(meta), path("tumor_fraction_estimation"), emit: results
    tuple val("${task.process}"), val('tumor_fraction_estimator'), eval('tumor_fraction_estimator --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args ?: ''
    def included_regions_opt = included_regions ? "--included-regions-bed '${included_regions}'" : "--included-regions-bed '/resources/GRCh38_GIAB_highconf.notinalldifficultregions.noY.bed'"
    def pon_opt             = panel_of_normals_tsv ? "--systematic-noise-tsv '${panel_of_normals_tsv}'" : ""
    def normal_sample_opt   = (meta_reporter_vcf.containsKey('normal_sample_name') && meta_reporter_vcf.normal_sample_name) ? "--normal-sample-name '${meta_reporter_vcf.normal_sample_name}'" : ""

    """
    mkdir -p 'tumor_fraction_estimation'

    tumor_fraction_estimator tumor-fraction \\
        --bam '${bam}' \\
        --query-vcf '${vcf}' \\
        --tumor-sample-name '${meta_reporter_vcf.tumor_sample_name}' \\
        ${normal_sample_opt} \\
        --reference '${reference}' \\
        --output-dir 'tumor_fraction_estimation' \\
        ${included_regions_opt} \\
        ${pon_opt} \\
        --threads ${task.cpus} \\
        ${args}
    """

    stub:
    """
    mkdir -p tumor_fraction_estimation/metrics
    touch "tumor_fraction_estimation/tumor_fraction_estimation_results.tsv"
    touch "tumor_fraction_estimation/metrics/per_site_counts.tsv"
    """
}
