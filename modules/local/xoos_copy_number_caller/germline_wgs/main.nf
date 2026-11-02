include { libraryTypeFromRunType } from '../../../../subworkflows/local/utils_xoos_workflow/main'

process XOOS_COPY_NUMBER_CALLER_GERMLINE_WGS {
    conda false
    container 'ghcr.io/roche-axelios/xoos/copy_number_caller:1.2.0'

    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index), path(small_variant_vcf), path(small_variant_vcf_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(mappability_bigwig)
    tuple val(meta4), path(cnv_seed_file)
    tuple val(meta5), path(blocklist_bed)

    output:
    tuple val(meta), path("copy_number_caller/${prefix}_germline_copy_number_callset.vcf.gz"), emit: vcf
    tuple val(meta), path("copy_number_caller/*"), emit: results_dir
    tuple val("${task.process}"), val('copy_number_caller'), eval("copy_number_caller --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def library_type = libraryTypeFromRunType(meta.run_type)
    def vcf_arg = small_variant_vcf ? "--vcf '${small_variant_vcf}'" : ''
    def mappability_arg = mappability_bigwig ? "--mappability-bigwig '${mappability_bigwig}'" : ''
    // HaplotypeCaller derives the VCF sample column from the BAM SM tag. CNC
    // must use that same value when a BAM's SM differs from the samplesheet ID.
    def sample_name = meta.bam_sample_name ?: meta.germline_sample_name ?: meta.sample_name ?: meta.id
    def seg_path = "copy_number_caller/${prefix}_germline_copy_number_callset.seg"
    def reheader_seg = sample_name != meta.id ? """
    # CNC uses --sample-name as the SEG sample identifier too. Normalize the
    # SEG sample column after the call so all published CNC artifacts use the
    # canonical samplesheet ID.
    if [ -f '${seg_path}' ]; then
        awk -v old='${sample_name}' -v new='${meta.id}' \\
            'BEGIN { FS = OFS = "\\t" } \$1 == old { \$1 = new } { print }' \\
            '${seg_path}' > '${seg_path}.tmp'
        mv '${seg_path}.tmp' '${seg_path}'
    fi
    """ : ''
    """
    mkdir -p "copy_number_caller"

    copy_number_caller germline-wgs \\
        --bam '${bam}' \\
        --reference '${reference}' \\
        ${mappability_arg} \\
        --seed-segments '${cnv_seed_file ?: '/resources/hg38_seed.seg'}' \\
        --blocklist-bed '${blocklist_bed ?: '/resources/hg38_decoy_intervals.bed'}' \\
        --sample-name '${sample_name}' \\
        --library-type '${library_type}' \\
        --output-dir 'copy_number_caller' \\
        --output-prefix '${prefix}' \\
        ${vcf_arg} \\
        --threads ${task.cpus} \\
        ${args}
    ${reheader_seg}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = meta.bam_sample_name ?: meta.germline_sample_name ?: meta.sample_name ?: meta.id
    def seg_path = "copy_number_caller/${prefix}_germline_copy_number_callset.seg"
    def reheader_seg = sample_name != meta.id ? """
    awk -v old='${sample_name}' -v new='${meta.id}' \\
        'BEGIN { FS = OFS = "\\t" } \$1 == old { \$1 = new } { print }' \\
        '${seg_path}' > '${seg_path}.tmp'
    mv '${seg_path}.tmp' '${seg_path}'
    """ : ''
    def metrics_dir = 'copy_number_caller/metrics'
    def tracks_dir = 'copy_number_caller/tracks'
    def supplemental_dir = 'copy_number_caller/supplemental'
    """
    mkdir -p '${metrics_dir}' '${tracks_dir}' '${supplemental_dir}'
    printf 'sample\\tchromosome\\tstart\\tend\\n${sample_name}\\tchr1\\t1\\t100\\n' > '${seg_path}'
    ${reheader_seg}
    gzip -n < /dev/null > 'copy_number_caller/${prefix}_germline_copy_number_callset.vcf.gz'
    touch 'copy_number_caller/${prefix}_germline_copy_number_callset.vcf.gz.tbi'
    touch '${metrics_dir}/${prefix}_metrics.tsv'
    touch '${tracks_dir}/${prefix}_igv_visualization.xml'
    touch '${tracks_dir}/${prefix}_log_ratios.bw'
    touch '${tracks_dir}/${prefix}_bafs.bw'
    touch '${supplemental_dir}/${prefix}_augmented_intervals.bed'
    touch '${supplemental_dir}/${prefix}_coverage.bed'
    touch '${supplemental_dir}/${prefix}_corrected_coverage.bed'
    touch '${supplemental_dir}/${prefix}_log_ratios.bed'
    touch '${supplemental_dir}/${prefix}_mapping_quality.bed'
    touch '${supplemental_dir}/${prefix}_bafs.bed'
    """
}
