include { XOOS_READ_COLLAPSER_RESCUE_SECONDARY } from '../../../modules/local/xoos_read_collapser/rescue_secondary/main'
include { XOOS_READ_COLLAPSER_MARK_DUPLICATES  } from '../../../modules/local/xoos_read_collapser/markdup/main'
include { XOOS_PARSE_READ_COLLAPSER_METRICS    } from '../../../modules/local/xoos_generate_metrics_report/parse_read_collapser_metrics/main'
include { XOOS_PARSE_RESCUE_SECONDARY_METRICS  } from '../../../modules/local/xoos_generate_metrics_report/parse_rescue_secondary_metrics/main'
include { SAMTOOLS_INDEX                       } from '../../../modules/nf-core/samtools/index/main'

workflow MARKDUP {
    take:
    bams // channel: [ val(meta), path(bam), path(bai) ]
    enable_secondary_rescue // val: true if this run needs secondary-alignment rescue

    main:
    // 1. Optionally rescue secondary alignments (for pangenome-aligned data).
    ch_rescue_secondary_report = channel.empty()
    ch_rescued = channel.empty()
    if (enable_secondary_rescue) {
        bams
            .branch { meta, _bam, _bai ->
                rescue: !meta.bam_secondary_rescued
                bypass: true
            }
            .set { ch_rescue_branch }

        ch_rescued_unindexed = XOOS_READ_COLLAPSER_RESCUE_SECONDARY(ch_rescue_branch.rescue).rescued_bam

        samtools_index_out = SAMTOOLS_INDEX(ch_rescued_unindexed)

        ch_rescued = ch_rescued_unindexed
            .join(samtools_index_out.index)
            .map { meta, bam, bai -> [meta + [bam_secondary_rescued: true], bam, bai] }
            .mix(ch_rescue_branch.bypass)

        ch_rescue_secondary_report = XOOS_PARSE_RESCUE_SECONDARY_METRICS(
            XOOS_READ_COLLAPSER_RESCUE_SECONDARY.out.metrics_dir
        ).report
    }

    // 2. Run read_collapser markdup. Every rescued BAM needs to run through markdup
    // even if it was previously marked.
    rc_out = XOOS_READ_COLLAPSER_MARK_DUPLICATES(enable_secondary_rescue ? ch_rescued : bams)

    // 3. Join the BAM and BAI back together into the final shape
    ch_markdup_bam = rc_out.markdup_bams
        .map { meta, bam, bai, _outdir -> [meta, bam, bai] }

    // 4. Parse the metrics (markdup + rescue_secondary if enabled)
    ch_metrics_report = XOOS_PARSE_READ_COLLAPSER_METRICS(rc_out.metrics_dir).report
        .mix(ch_rescue_secondary_report)

    emit:
    markdup_bam                = ch_markdup_bam  // channel: [ val(meta), path(bam), path(bai) ]
    report                     = ch_metrics_report // channel: [ val(meta), path(json_report) ]
    read_collapser_metrics_dir = rc_out.metrics_dir // channel: [ val(meta), path(metrics_dir) ]
    versions                   = Channel.empty() // versions are collected via topic channels
}
