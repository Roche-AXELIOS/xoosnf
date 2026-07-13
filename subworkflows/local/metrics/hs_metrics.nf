include { PICARD_COLLECTHSMETRICS         } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/picard/createsequencedictionary/main'
include { XOOS_PARSE_HS_METRICS           } from '../../../modules/local/xoos_generate_metrics_report/parse_hs_metrics/main'

workflow METRICS_HS {
    take:
    name
    bams // channel: [ val(meta), path(bam), path(bai) ]
    ch_picard_reference // channel: [ val(meta), path(fasta), path(fai) ]
    ch_capture_bed // channel: [ val(meta), path(bed) ]
    ch_target_bed // channel: [ val(meta), path(bed) ]

    main:

    // Create sequence dictionary for reference
    def ref_meta = [id: 'reference']
    def ch_reference_dict = PICARD_CREATESEQUENCEDICTIONARY(ch_picard_reference.map { meta, fasta, _fai -> [meta, fasta] }).reference_dict
    // channel: [ val(meta), path(fasta) ]

    // Prepare HS metrics inputs with BED files
    hs_inputs = bams
        .combine(ch_capture_bed)
        .combine(ch_target_bed)
        .map { meta, bam, bai, _capture_bed_meta, capture_bed, _target_bed_meta, target_bed ->
            // Use per-sample BED file specified in metadata if available, otherwise use the one from the channel
            def capture_bed_final = meta.containsKey('capture_bed') ? file(meta.capture_bed) : capture_bed
            def target_bed_final = meta.containsKey('target_bed') ? file(meta.target_bed) : target_bed
            if (!capture_bed_final || !target_bed_final) {
                error "Sample ${meta.id}: HS metrics enabled but capture and/or target BED is missing. Provide capture_bed and target_bed, or set enable_hs_metrics to false."
            }
            return [meta, bam, bai, capture_bed_final, target_bed_final]
        }

    // Collect HS metrics and generate report
    ch_hs_metrics = PICARD_COLLECTHSMETRICS(
        hs_inputs,
        ch_picard_reference.map { meta, fasta, _fai -> [meta, fasta] },
        ch_picard_reference.map { meta, _fasta, fai -> [meta, fai] },
        ch_reference_dict,
        channel.value([ref_meta, []]),
    ).metrics
    hs_metrics_report = XOOS_PARSE_HS_METRICS(
        ch_hs_metrics.map { meta, metrics -> [meta, name, metrics] }
    ).report

    emit:
    report        = hs_metrics_report // channel: [ val(meta), path(report) ]
    metrics_files = ch_hs_metrics // channel: [ val(meta), path(metrics_dir) ]
    versions      = Channel.empty() // versions are collected via topic channels
}
