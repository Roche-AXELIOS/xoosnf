include { XOOS_READ_COLLAPSER_RESCUE_SECONDARY } from '../../../modules/local/xoos_read_collapser/rescue_secondary/main'
include { XOOS_READ_COLLAPSER_MARK_DUPLICATES  } from '../../../modules/local/xoos_read_collapser/markdup/main'
include { XOOS_PARSE_READ_COLLAPSER_METRICS    } from '../../../modules/local/xoos_generate_metrics_report/parse_read_collapser_metrics/main'
include { XOOS_PARSE_RESCUE_SECONDARY_METRICS  } from '../../../modules/local/xoos_generate_metrics_report/parse_rescue_secondary_metrics/main'
include { SAMTOOLS_INDEX                       } from '../../../modules/nf-core/samtools/index/main'

workflow MARKDUP {
    take:
    bams // channel: [ val(meta), path(bam), path(bai) ]
    rescue_policy // val: true if this run's configuration wants secondary-alignment rescue

    main:
    // 1. Optionally rescue secondary alignments (for pangenome-aligned data).
    //    The decision is per-BAM: rescue runs only when the run wants it AND the
    //    BAM has not already been rescued. Secondary rescue is not idempotent —
    //    re-running it on an already-rescued BAM (e.g. a re-fed pipeline output)
    //    would corrupt the alignment record set, so already-rescued BAMs bypass
    //    the step while still proceeding to markdup.
    ch_rescue_secondary_report = channel.empty()
    if (rescue_policy) {
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
            // Stamp the rescued state. This flag is preserved through later steps
            // and prevents a second rescue if the BAM ever re-enters the pipeline.
            .map { meta, bam, bai -> [meta + [bam_secondary_rescued: true], bam, bai] }
            .mix(ch_rescue_branch.bypass)

        ch_rescue_secondary_report = XOOS_PARSE_RESCUE_SECONDARY_METRICS(
            XOOS_READ_COLLAPSER_RESCUE_SECONDARY.out.metrics_dir
        ).report
    } else {
        ch_rescued = bams
    }

    // 2. Run read_collapser markdup
    rc_out = XOOS_READ_COLLAPSER_MARK_DUPLICATES(ch_rescued)

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
