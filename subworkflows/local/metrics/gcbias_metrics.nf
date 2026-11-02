include { XOOS_ALIGNMENT_METRICS_GC_BIAS       } from '../../../modules/local/xoos_alignment_metrics/gcbias/main'
include { XOOS_PARSE_ALIGNMENT_METRICS_GC_BIAS } from '../../../modules/local/xoos_generate_metrics_report/parse_alignment_metrics_gc_bias/main'
include { XOOS_GC_BIAS_METRICS_TO_MULTIQC      } from '../../../modules/local/xoos_metrics_to_multiqc/gc_bias/main'

workflow METRICS_GC_BIAS {
    take:
    name
    bams // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    // Preserve the metrics stage on the GC-bias channel. This mirrors the
    // alignment-metrics path and lets the MultiQC conversion keep results from
    // the primary and downsampled alignments separate.
    ch_gc_bias_metrics = XOOS_ALIGNMENT_METRICS_GC_BIAS(
        bams,
        ch_reference,
    ).gc_bias_metrics.map { meta, metrics -> [meta + [report_prefix: name], metrics] }
    gc_bias_metrics_report = XOOS_PARSE_ALIGNMENT_METRICS_GC_BIAS(
        ch_gc_bias_metrics.map { meta, metrics -> [meta, name, metrics] }
    ).report
    ch_multiqc_files = XOOS_GC_BIAS_METRICS_TO_MULTIQC(
        ch_gc_bias_metrics
    ).multiqc_report_dir

    emit:
    report        = gc_bias_metrics_report // channel: [ val(meta), path(report) ]
    metrics_files = ch_gc_bias_metrics // channel: [ val(meta), path(metrics_dir) ]
    multiqc_files = ch_multiqc_files // channel: [ val(meta), path(multiqc_input_dir) ]
    versions      = Channel.empty() // versions are collected via topic channels
}
