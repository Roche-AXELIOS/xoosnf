include { XOOS_ALIGNMENT_METRICS_ALL as ALIGNMENT_METRICS_AUTOSOMES               } from '../../../modules/local/xoos_alignment_metrics/all/main'
include { XOOS_ALIGNMENT_METRICS_ALL as ALIGNMENT_METRICS_TARGETS                 } from '../../../modules/local/xoos_alignment_metrics/all/main'
include { XOOS_ALIGNMENT_METRICS_ALL as ALIGNMENT_METRICS_HIGH_CONFIDENCE         } from '../../../modules/local/xoos_alignment_metrics/all/main'
include { XOOS_PARSE_ALIGNMENT_METRICS as PARSE_ALIGNMENT_METRICS_AUTOSOMES       } from '../../../modules/local/xoos_generate_metrics_report/parse_alignment_metrics/main'
include { XOOS_PARSE_ALIGNMENT_METRICS as PARSE_ALIGNMENT_METRICS_TARGETS         } from '../../../modules/local/xoos_generate_metrics_report/parse_alignment_metrics/main'
include { XOOS_PARSE_ALIGNMENT_METRICS as PARSE_ALIGNMENT_METRICS_HIGH_CONFIDENCE } from '../../../modules/local/xoos_generate_metrics_report/parse_alignment_metrics/main'
include { XOOS_ALIGNMENT_METRICS_TO_MULTIQC                                       } from '../../../modules/local/xoos_metrics_to_multiqc/alignment_metrics/main'

workflow METRICS_ALIGNMENT {
    take:
    prefix // val: prefix for the report (e.g. "pre_dedup", "dedup", "dedup_subsample")
    ch_bam // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference // channel: [ val(meta), path(fasta), path(fai) ]
    ch_autosomes_bed // channel: [ val(meta), path(bed) ]
    ch_target_bed // channel: [ val(meta), path(bed) ]
    ch_high_confidence_bed_dir // channel: [ val(meta), path(bundle) ]
    ch_fallback_high_confidence_bed // channel: [ val(meta), path(bed) ]
    enable_alignment_metrics_autosomes
    enable_alignment_metrics_targets
    enable_alignment_metrics_high_confidence

    main:

    ch_report = channel.empty()

    // Add prefix to meta so tools can use it for naming outputs
    ch_bam = ch_bam.map { meta, bam, bai -> [meta + [report_prefix: prefix], bam, bai] }

    // Autosomes metrics
    ch_autosomes_metrics = channel.empty()
    if (enable_alignment_metrics_autosomes) {

        alignment_metrics_autosomes_out = ALIGNMENT_METRICS_AUTOSOMES(
            ch_bam,
            ch_reference,
            ch_autosomes_bed,
        )
        ch_autosomes_metrics = alignment_metrics_autosomes_out.read_metrics
            .join(alignment_metrics_autosomes_out.coverage_metrics)
            .join(alignment_metrics_autosomes_out.accuracy_metrics)
            .map { meta, read_metrics, coverage, accuracy -> [meta + [metrics_region: 'autosomes'], read_metrics, coverage, accuracy] }

        ch_autosomes_report = PARSE_ALIGNMENT_METRICS_AUTOSOMES(
            ch_autosomes_metrics,
            prefix + '_autosomes',
        ).report
        ch_report = ch_report.concat(ch_autosomes_report)
    }

    // Target metrics
    ch_target_metrics = channel.empty()
    if (enable_alignment_metrics_targets) {
        ch_per_sample_target_bed = ch_bam
            .combine(ch_target_bed)
            .map { meta, _bam, _bai, _target_bed_meta, target_bed ->
                // Use per-sample BED file specified in metadata if available, otherwise use the one from the channel
                if (meta.containsKey('target_bed')) {
                    return [meta, file(meta.target_bed)]
                }
                else {
                    if (!target_bed) {
                        log.warn "Sample ${meta.id} is missing target BED file. Skipping target metrics for this sample."
                        return [meta, null]
                    }
                    return [meta, target_bed]
                }
            }
            .filter { _meta, target_bed -> target_bed != null }
        alignment_metrics_targets_out = ALIGNMENT_METRICS_TARGETS(
            ch_bam,
            ch_reference,
            ch_per_sample_target_bed,
        )
        ch_target_metrics = alignment_metrics_targets_out.read_metrics
            .join(alignment_metrics_targets_out.coverage_metrics)
            .join(alignment_metrics_targets_out.accuracy_metrics)
            .map { meta, read_metrics, coverage, accuracy -> [meta + [metrics_region: 'targets'], read_metrics, coverage, accuracy] }

        ch_target_report = PARSE_ALIGNMENT_METRICS_TARGETS(
            ch_target_metrics,
            prefix + '_targets',
        ).report
        ch_report = ch_report.concat(ch_target_report)
    }

    ch_high_confidence_metrics = channel.empty()
    if (enable_alignment_metrics_high_confidence) {
        // Match each BAM file with its corresponding high confidence BED file based on sample name
        ch_high_confidence_bed = ch_bam
            .combine(ch_high_confidence_bed_dir.map { _meta, bundle -> [bundle] })
            .combine(ch_fallback_high_confidence_bed.map { _meta, bed -> [bed] })
            .map { meta, _bam, _bai, high_confidence_bed_dir, fallback_bed ->
                // For each sample, find the corresponding high confidence BED file
                def bed_matching_sample_name = file(high_confidence_bed_dir)
                    .listFiles()
                    .find { file -> file.getName().contains(meta.id.split('_')[0]) && file.getName().endsWith('.bed') }
                if (bed_matching_sample_name) {
                    return [meta + [is_fallback: false], bed_matching_sample_name]
                }
                else {
                    // If no matching BED file is found, use the fallback BED file
                    return [meta + [is_fallback: true], file(fallback_bed)]
                }
            }
        alignment_metrics_high_confidence_out = ALIGNMENT_METRICS_HIGH_CONFIDENCE(
            ch_bam,
            ch_reference,
            ch_high_confidence_bed,
        )
        ch_high_confidence_metrics = alignment_metrics_high_confidence_out.read_metrics
            .join(alignment_metrics_high_confidence_out.coverage_metrics)
            .join(alignment_metrics_high_confidence_out.accuracy_metrics)
            .map { meta, read_metrics, coverage, accuracy -> [meta + [metrics_region: 'high_confidence'], read_metrics, coverage, accuracy] }

        ch_high_confidence_report = PARSE_ALIGNMENT_METRICS_HIGH_CONFIDENCE(
            ch_high_confidence_metrics,
            prefix + '_high_confidence',
        ).report
        ch_report = ch_report.concat(ch_high_confidence_report)
    }

    ch_alignment_metrics = ch_autosomes_metrics
        .concat(ch_high_confidence_metrics)
        .concat(ch_target_metrics)

    // Feed every enabled region into the MultiQC conversion so the report
    // includes all regions (autosomes / high_confidence / targets). Each item
    // carries meta.metrics_region and meta.report_prefix (the pipeline stage),
    // which the conversion turns into a selectable region/stage tab.
    ch_multiqc_conversion_input_files = ch_alignment_metrics
    ch_multiqc_files = XOOS_ALIGNMENT_METRICS_TO_MULTIQC(ch_multiqc_conversion_input_files).multiqc_report_dir

    emit:
    report        = ch_report // channel: [ val(meta), path(report) ]
    metrics_files = ch_alignment_metrics // channel: [ val(meta), path(read_metrics), path(coverage_metrics), path(accuracy_metrics) ]
    multiqc_files = ch_multiqc_files // channel: [ val(run_meta), path(multiqc_input_dir) ]
    versions      = Channel.empty() // versions are collected via topic channels
}
