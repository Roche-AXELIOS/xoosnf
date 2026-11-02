include { ALIGN                             } from '../align/main'

include { XOOS_READ_COLLAPSER_CONSENSUS     } from '../../../modules/local/xoos_read_collapser/consensus/main'
include { XOOS_PARSE_READ_COLLAPSER_METRICS } from '../../../modules/local/xoos_generate_metrics_report/parse_read_collapser_metrics/main'
workflow CONSENSUS_AND_REALIGN {
    take:
    ch_bam // channel: [ val(meta), path(bam), path(bai) ]
    ch_bed // channel: [ val(meta), path(bed) ]
    aligner // val: 'bwa', 'minimap2', 'vg_giraffe', 'pb_fq2bam', 'pb_minimap2', or 'pb_giraffe'
    ch_fasta // channel: [ val(meta), path(fasta), path(fai) ]
    ch_bwa_idx // channel: [ val(meta), path(idx_dir) ]
    ch_giraffe_idx // channel: [ val(meta), path(idx_dir), val(idx_name) ]
    ch_ref_paths // channel: [ path(dict) ] — linearized reference .dict for pbrun giraffe --ref-paths

    main:
    consensus_fastq = XOOS_READ_COLLAPSER_CONSENSUS(
        ch_bam.combine(ch_bed).map { meta, bam, bai, _meta2, bed -> [meta, bam, bai, bed] },
    ).fastq
    align_out = ALIGN(
        aligner,
        consensus_fastq,
        ch_fasta,
        ch_bwa_idx,
        ch_giraffe_idx,
        ch_ref_paths,
    )
    consensus_bam = align_out.bams

    ch_read_collapser_metrics_json_report = XOOS_PARSE_READ_COLLAPSER_METRICS(XOOS_READ_COLLAPSER_CONSENSUS.out.metrics_dir).report

    emit:
    bam                        = consensus_bam // channel: [ val(meta), path(bam), path(bai) ]
    giraffe_graphs             = align_out.giraffe_graphs // channel: [ val(meta), path(gbz), path(zipcodes), path(min), path(dist) ]
    report                     = ch_read_collapser_metrics_json_report // channel: [ val(meta), path(json_report) ]
    read_collapser_metrics_dir = XOOS_READ_COLLAPSER_CONSENSUS.out.metrics_dir // channel: [ val(meta), path(metrics_dir) ]
    versions                   = Channel.empty() // versions are collected via topic channels
}
