include { XOOS_TUMOR_FRACTION_ESTIMATOR_TUMOR_FRACTION } from '../../../modules/local/xoos_tumor_fraction_estimator/tumor_fraction/main'
include { XOOS_TFE_METRICS_TO_MULTIQC                  } from '../../../modules/local/xoos_metrics_to_multiqc/tfe_metrics/main'

workflow TUMOR_FRACTION_ESTIMATION {

    take:
    ch_bams                   // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference              // channel: [ val(meta_ref), path(fasta), path(fai) ]
    ch_tfe_included_regions   // channel: [ val(meta_bed), path(bed) ]
    ch_tfe_panel_of_normals     // channel: [ val(meta_pon), path(panel_of_normals_tsv), path(panel_of_normals_tbi) ]

    main:

    ch_custom_multiqc_files = channel.empty()

    // TFE input comes from reporter_vcf on the samplesheet.
    // The reporter_vcf may originate from a separate somatic TN run or be
    // provided directly for cfDNA workflows.
    // Samples without reporter_vcf are excluded — TFE produces no output for them.
    ch_tfe_samples = ch_bams
        .filter { meta, _bam, _bai -> meta.containsKey('reporter_vcf') && meta.reporter_vcf }

    ch_tfe_bam_input = ch_tfe_samples
        .map { meta, bam, bai -> [ meta, bam, bai ] }

    ch_reporter_vcf = ch_tfe_samples
        .map { meta, _bam, _bai ->
            def meta_reporter_vcf = [ tumor_sample_name: meta.reporter_vcf_tumor_sample ]
            if (meta.containsKey('reporter_vcf_normal_sample') && meta.reporter_vcf_normal_sample) {
                meta_reporter_vcf = meta_reporter_vcf + [ normal_sample_name: meta.reporter_vcf_normal_sample ]
            }
            [ meta_reporter_vcf, file(meta.reporter_vcf), [] ]
        }

    XOOS_TUMOR_FRACTION_ESTIMATOR_TUMOR_FRACTION(
        ch_tfe_bam_input,
        ch_reporter_vcf,
        ch_reference,
        ch_tfe_included_regions,
        ch_tfe_panel_of_normals
    )

    XOOS_TFE_METRICS_TO_MULTIQC(XOOS_TUMOR_FRACTION_ESTIMATOR_TUMOR_FRACTION.out.results)
    ch_custom_multiqc_files = ch_custom_multiqc_files.mix(XOOS_TFE_METRICS_TO_MULTIQC.out.multiqc_report_dir)

    emit:
    results              = XOOS_TUMOR_FRACTION_ESTIMATOR_TUMOR_FRACTION.out.results
    custom_multiqc_files = ch_custom_multiqc_files
}
