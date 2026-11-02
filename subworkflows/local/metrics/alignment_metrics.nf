include { XOOS_MERGE_AND_LABEL_BED                                                     } from '../../../modules/local/xoos_alignment_metrics/merge_and_label_bed/main'
include { XOOS_ALIGNMENT_METRICS_ALL_STRATIFIED as ALIGNMENT_METRICS_STRATIFIED        } from '../../../modules/local/xoos_alignment_metrics/all_stratified/main'
include { XOOS_ALIGNMENT_METRICS_ALL as ALIGNMENT_METRICS_HIGH_CONFIDENCE              } from '../../../modules/local/xoos_alignment_metrics/all/main'
include { XOOS_PARSE_ALIGNMENT_METRICS as PARSE_ALIGNMENT_METRICS                      } from '../../../modules/local/xoos_generate_metrics_report/parse_alignment_metrics/main'
include { XOOS_PARSE_ALIGNMENT_METRICS as PARSE_ALIGNMENT_METRICS_HIGH_CONFIDENCE      } from '../../../modules/local/xoos_generate_metrics_report/parse_alignment_metrics/main'
include { XOOS_ALIGNMENT_METRICS_TO_MULTIQC                                            } from '../../../modules/local/xoos_metrics_to_multiqc/alignment_metrics/main'

workflow METRICS_ALIGNMENT {
    take:
    prefix // val: stage prefix for the report (e.g. "pre_dedup", "dedup", "dedup_subsample", "downsampled_metrics")
    aligner // val: aligner name inserted into the metric key (e.g. "bwa", "pb_fq2bam")
    ch_bam // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference // channel: [ val(meta), path(fasta), path(fai) ]
    ch_autosomes_bed // channel: [ val(meta), path(bed) ]
    ch_target_bed // channel: [ val(meta), path(bed) ]
    ch_high_confidence_bed_dir // channel: [ val(meta), path(bundle) ]
    ch_fallback_high_confidence_bed // channel: [ val(meta), path(bed) ]
    enable_alignment_metrics_autosomes
    enable_alignment_metrics_targets
    enable_alignment_metrics_high_confidence
    additional_bed_files // list: [ [region_name: '<label>', bed: '<path>'], ... ]

    main:

    ch_report = channel.empty()

    ch_bam = ch_bam.map { meta, bam, bai -> [meta + [report_prefix: prefix], bam, bai] }

    // All non-HC region BEDs are merged into a single labeled BED per sample,
    // then run through one ALIGNMENT_METRICS_STRATIFIED call. The output is
    // fanned out per-label for parsing and MultiQC.

    def extra_labels = additional_bed_files.collect { it -> it.region_name }
    def extra_beds = additional_bed_files.collect { it -> file(it.bed, checkIfExists: true) }

    ch_stratified_metrics = channel.empty()

    if (enable_alignment_metrics_autosomes || enable_alignment_metrics_targets || additional_bed_files) {

        ch_merge_input = ch_bam
            .combine(ch_autosomes_bed.map { _m, bed -> [bed] })
            .combine(ch_target_bed.map { _m, bed -> [bed] })
            .map { meta, _bam, _bai, autosomes_bed, target_bed ->
                def labels = extra_labels.collect()
                def beds = extra_beds.collect()

                if (enable_alignment_metrics_autosomes && autosomes_bed) {
                    labels << 'autosomes'
                    beds << autosomes_bed
                }

                if (enable_alignment_metrics_targets) {
                    def per_sample_target = meta.containsKey('target_bed') ? file(meta.target_bed) : target_bed
                    if (per_sample_target) {
                        labels << 'targets'
                        beds << per_sample_target
                    } else {
                        log.warn "Sample ${meta.id}: no target BED available, skipping targets region."
                    }
                }
                return [meta + [region_labels: labels], labels, beds]
            }
            .filter { _m, labels, _beds -> !labels.isEmpty() }

        ch_labeled_bed = XOOS_MERGE_AND_LABEL_BED(ch_merge_input).labeled_bed

        // Inject region_labels into the BAM meta so the stratified process output
        // carries the list of region labels for each sample for downstream fan-out.
        ch_region_labels = ch_labeled_bed.map { meta, _bed -> [meta.id, meta.region_labels] }
        ch_bam_with_labels = ch_bam
            .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
            .combine(ch_region_labels, by: 0)
            .map { _id, meta, bam, bai, region_labels -> [meta + [region_labels: region_labels], bam, bai] }

        alignment_metrics_out = ALIGNMENT_METRICS_STRATIFIED(
            ch_bam_with_labels,
            ch_reference,
            ch_labeled_bed,
        )
        // Fan out the stratified metrics by region labels
        // One entry per region label for each sample
        ch_stratified_metrics = alignment_metrics_out.metrics_dir.flatMap { meta, dir ->
            meta.region_labels.collect { label ->
                [
                    meta + [metrics_region: label],
                    dir.resolve(label).resolve('read_metrics'),
                    dir.resolve(label).resolve('coverage'),
                    dir.resolve(label).resolve('accuracy'),
                ]
            }
        }

        ch_parse_input = ch_stratified_metrics.multiMap { meta, read_metrics, coverage, accuracy ->
            metrics: [meta, read_metrics, coverage, accuracy]
            source: "${prefix}_${meta.metrics_region}"
            region: "${prefix}/${aligner}/${meta.metrics_region}"
        }
        ch_report = ch_report.concat(
            PARSE_ALIGNMENT_METRICS(
                ch_parse_input.metrics,
                ch_parse_input.region,
                ch_parse_input.source,
            ).report
        )
    }

    // Run high-confidence metrics for each sample using the appropriate BED files.
    // The high-confidence metrics are generated separately from the stratified metrics
    // because they use a different `--max-alt-allele-fraction`.

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
                    .find { f -> f.getName().contains(meta.id.split('_')[0]) && f.getName().endsWith('.bed') }
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
            "${prefix}/${aligner}/high_confidence",
            "${prefix}_high_confidence",
        ).report
        ch_report = ch_report.mix(ch_high_confidence_report)
    }

    // Combine all metrics for MultiQC
    ch_alignment_metrics = ch_stratified_metrics
        .mix(ch_high_confidence_metrics)

    // Feed every enabled region into the MultiQC conversion so the report
    // includes all regions (autosomes, high_confidence, targets, etc). Each item
    // carries meta.metrics_region and meta.report_prefix (the pipeline stage),
    // which the conversion turns into a selectable region/stage tab.
    ch_multiqc_files = XOOS_ALIGNMENT_METRICS_TO_MULTIQC(ch_alignment_metrics).multiqc_report_dir

    emit:
    report        = ch_report // channel: [ val(meta), path(report) ]
    metrics_files = ch_alignment_metrics // channel: [ val(meta), path(read_metrics), path(coverage_metrics), path(accuracy_metrics) ]
    multiqc_files = ch_multiqc_files // channel: [ val(run_meta), path(multiqc_input_dir) ]
    versions      = channel.empty() // versions are collected via topic channels
}
