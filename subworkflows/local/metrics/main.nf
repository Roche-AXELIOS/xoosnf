include { METRICS_ALIGNMENT } from './alignment_metrics'
include { METRICS_GC_BIAS   } from './gcbias_metrics'
include { METRICS_HS        } from './hs_metrics'

workflow METRICS {
    take:
    name
    ch_bam
    ch_autosomes_bed
    ch_target_bed
    ch_capture_bed
    ch_high_confidence_bed_dir
    ch_fallback_high_confidence_bed
    ch_reference_for_metrics
    enable_alignment_metrics_autosomes
    enable_alignment_metrics_targets
    enable_alignment_metrics_high_confidence
    enable_gc_bias_metrics
    enable_hs_metrics
    aligner
    additional_bed_files // list: [ [region_name, bed] ]

    main:
    ch_versions = Channel.empty()

    ch_multiqc_files = channel.empty()
    ch_custom_multiqc_files = channel.empty()

    ch_alignment_metrics_report = METRICS_ALIGNMENT(
        name,
        aligner,
        ch_bam,
        ch_reference_for_metrics,
        ch_autosomes_bed,
        ch_target_bed,
        ch_high_confidence_bed_dir,
        ch_fallback_high_confidence_bed,
        enable_alignment_metrics_autosomes,
        enable_alignment_metrics_targets,
        enable_alignment_metrics_high_confidence,
        additional_bed_files,
    ).report

    ch_alignment_metrics_files = METRICS_ALIGNMENT.out.metrics_files
    ch_custom_multiqc_files = ch_custom_multiqc_files.mix(METRICS_ALIGNMENT.out.multiqc_files)
    ch_versions = ch_versions.mix(METRICS_ALIGNMENT.out.versions)

    // Run GC bias metrics subworkflow
    ch_gc_bias_metrics_report = channel.empty()
    ch_gc_bias_metrics_files = channel.empty()
    if (enable_gc_bias_metrics) {
        ch_gc_bias_metrics_report = METRICS_GC_BIAS(
            name,
            ch_bam,
            ch_reference_for_metrics,
        ).report
        ch_gc_bias_metrics_files = METRICS_GC_BIAS.out.metrics_files
        ch_custom_multiqc_files = ch_custom_multiqc_files.mix(METRICS_GC_BIAS.out.multiqc_files)
        ch_versions = ch_versions.mix(METRICS_GC_BIAS.out.versions)
    }

    // Run Hybrid Selection (HS) metrics subworkflow
    ch_hs_metrics_report = channel.empty()
    ch_hs_metrics_files = channel.empty()
    if (enable_hs_metrics) {
        ch_hs_metrics_report = METRICS_HS(
            name,
            ch_bam,
            ch_reference_for_metrics,
            ch_capture_bed,
            ch_target_bed,
        ).report
        ch_hs_metrics_files = METRICS_HS.out.metrics_files
            .mix(METRICS_HS.out.per_target_coverage)
            .mix(METRICS_HS.out.per_base_coverage)
        ch_multiqc_files = ch_multiqc_files.mix(METRICS_HS.out.metrics_files)
        ch_versions = ch_versions.mix(METRICS_HS.out.versions)
    }

    metrics_report = ch_alignment_metrics_report
        .concat(ch_gc_bias_metrics_report)
        .concat(ch_hs_metrics_report)

    emit:
    report                  = metrics_report
    alignment_metrics_files = ch_alignment_metrics_files // channel: [ val(meta), path(read_metrics), path(coverage), path(accuracy) ]
    gc_bias_metrics_files   = ch_gc_bias_metrics_files // channel: [ val(meta), path(gc_bias) ]
    hs_metrics_files        = ch_hs_metrics_files // channel: [ val(meta), path(*_metrics) ]
    custom_multiqc_files    = ch_custom_multiqc_files
    multiqc_files           = ch_multiqc_files
    versions                = ch_versions
}
