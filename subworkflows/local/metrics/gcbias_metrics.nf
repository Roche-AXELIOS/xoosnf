include { PICARD_COLLECT_GC_BIAS_METRICS } from '../../../modules/local/picard/collectgcbiasmetrics/main'
include { XOOS_PARSE_GC_BIAS_METRICS     } from '../../../modules/local/xoos_generate_metrics_report/parse_gc_bias/main'

workflow METRICS_GC_BIAS {
    take:
    name
    bams // channel: [ val(meta), path(bam), path(bai) ]
    ch_picard_reference // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    ch_gc_bias_metrics = PICARD_COLLECT_GC_BIAS_METRICS(
        bams,
        ch_picard_reference,
    ).gc_bias
    gc_bias_metrics_report = XOOS_PARSE_GC_BIAS_METRICS(
        ch_gc_bias_metrics.map { meta, metrics -> [meta, name, metrics] }
    ).report

    emit:
    report        = gc_bias_metrics_report // channel: [ val(meta), path(report) ]
    metrics_files = ch_gc_bias_metrics // channel: [ val(meta), path(metrics_dir) ]
    versions      = Channel.empty() // versions are collected via topic channels
}
