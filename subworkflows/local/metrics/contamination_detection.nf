include { XOOS_TUMOR_FRACTION_ESTIMATOR_CONTAMINATION_DETECTION } from '../../../modules/local/xoos_tumor_fraction_estimator/contamination_detection/main'
include { XOOS_CONTAMINATION_METRICS_TO_MULTIQC                 } from '../../../modules/local/xoos_metrics_to_multiqc/contamination_metrics/main'

workflow CONTAMINATION_DETECTION {

    take:
    ch_bams                          // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference                     // channel: [ val(meta_ref), path(fasta), path(fai) ]
    ch_contamination_vcf             // channel: [ val(meta_vcf), path(vcf) ]
    ch_contamination_included_regions // channel: [ val(meta_bed), path(bed) ]

    main:

    ch_custom_multiqc_files = channel.empty()

    // Extract the VCF and BED paths from their meta-wrapped channels.
    // When no external file is provided the path will be empty and the
    // binary uses its bundled default.
    ch_vcf_path = ch_contamination_vcf.map { _meta, vcf -> vcf ?: [] }
    ch_bed_path = ch_contamination_included_regions.map { _meta, bed -> bed ?: [] }

    XOOS_TUMOR_FRACTION_ESTIMATOR_CONTAMINATION_DETECTION(
        ch_bams,
        ch_reference,
        ch_vcf_path,
        ch_bed_path,
    )

    XOOS_CONTAMINATION_METRICS_TO_MULTIQC(XOOS_TUMOR_FRACTION_ESTIMATOR_CONTAMINATION_DETECTION.out.results)
    ch_custom_multiqc_files = ch_custom_multiqc_files.mix(XOOS_CONTAMINATION_METRICS_TO_MULTIQC.out.multiqc_report_dir)

    emit:
    results              = XOOS_TUMOR_FRACTION_ESTIMATOR_CONTAMINATION_DETECTION.out.results
    custom_multiqc_files = ch_custom_multiqc_files
}
