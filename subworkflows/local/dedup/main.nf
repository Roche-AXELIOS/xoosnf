include { MARKDUP               } from './markdup'
include { CONSENSUS_AND_REALIGN } from './consensus'

workflow DEDUP {
    take:
    dedup_strategy // val: 'markdup' or 'consensus'
    aligner // val: 'bwa', 'minimap2', 'vg_giraffe', 'pb_fq2bam', 'pb_minimap2', or 'pb_giraffe'
    enable_read_collapser_dedup // val: true or false — when true, the required dedup state is xoos_markdup regardless of aligner
    rescue_policy // val: true if this run's configuration wants secondary-alignment rescue (actual per-BAM decision is made in MARKDUP)
    ch_bam // channel: [ val(meta), path(bam), path(bai) ]
    ch_bed // channel: [ val(meta), path(bed) ]
    ch_fasta // channel: [ val(meta), path(fasta), path(fai) ]
    ch_bwa_idx // channel: [ val(meta), path(idx_dir) ]
    ch_giraffe_idx // channel: [ val(meta), path(idx_dir), val(idx_name) ]
    ch_ref_paths // channel: [ path(dict) ] — linearized reference .dict for pbrun giraffe --ref-paths

    main:
    // Collect versions (also published via topic channels)
    ch_versions = Channel.empty()

    dedup_bam = channel.empty()
    dedup_report = channel.empty()
    dedup_read_collapser_metrics_dir = channel.empty()
    if (dedup_strategy == 'markdup') {
        // The dedup states that are acceptable on an incoming BAM for this run,
        // i.e. states for which running MARKDUP again would be redundant. When
        // the read collapser dedup is enabled only its own marking (xoos_markdup)
        // is acceptable. For Parabricks profiles the aligner marks duplicates
        // (pb_markdup), a re-fed pipeline output (xoos_markdup) is also fully
        // marked, and external AXELIOS/SAS marking (sas_markdup) is trusted as a
        // valid pre-marked state — MARKDUP would only re-produce xoos_markdup —
        // so all three are accepted. When a run additionally requires secondary
        // rescue (rescue_policy), the per-BAM branch below still re-marks these
        // BAMs because rescue re-marks duplicates as part of the same step. Any
        // other profile requires xoos_markdup.
        def acceptable_dedup = enable_read_collapser_dedup ? ['xoos_markdup']
            : (aligner in ['pb_fq2bam', 'pb_minimap2', 'pb_giraffe']) ? ['sas_markdup', 'pb_markdup', 'xoos_markdup']
            : ['xoos_markdup']

        // Decide per-BAM whether the markdup subworkflow needs to run. It is
        // skipped only when the BAM is already in an acceptable dedup state AND
        // does not still need secondary rescue. Skipping a BAM that already
        // satisfies both avoids redundant markdup and, critically, prevents a
        // second (non-idempotent) secondary rescue on re-fed pipeline outputs.
        ch_bam
            .branch { meta, _bam, _bai ->
                def needs_rescue = rescue_policy && !meta.bam_secondary_rescued
                skip: meta.bam_dedup in acceptable_dedup && !needs_rescue
                run: true
            }
            .set { ch_markdup_branch }

        markdup_out = MARKDUP(ch_markdup_branch.run, rescue_policy)
        // Stamp the dedup state on freshly marked BAMs. bam_secondary_rescued is
        // updated inside MARKDUP (it owns the rescue step), so only bam_dedup is
        // set here. The rescued flag is preserved by the merge.
        dedup_bam = markdup_out.markdup_bam
            .map { meta, bam, bai -> [meta + [bam_dedup: 'xoos_markdup'], bam, bai] }
            .mix(ch_markdup_branch.skip)
        dedup_report = markdup_out.report
        dedup_read_collapser_metrics_dir = markdup_out.read_collapser_metrics_dir
        ch_versions = ch_versions.mix(markdup_out.versions)
    }
    else if (dedup_strategy == 'consensus') {
        // Secondary-alignment rescue is only wired into the markdup path. If the
        // run configuration requests rescue while using the consensus strategy,
        // it will not be applied — warn so this silent no-op is visible.
        if (rescue_policy) {
            log.warn("rescue_policy is enabled but dedup_strategy is 'consensus'; secondary-alignment rescue is only applied in the markdup path and will be skipped.")
        }
        consensus_out = CONSENSUS_AND_REALIGN(
            ch_bam,
            ch_bed,
            aligner,
            ch_fasta,
            ch_bwa_idx,
            ch_giraffe_idx,
            ch_ref_paths,
        )
        // Consensus re-aligns and collapses; stamp the resulting dedup state.
        // bam_secondary_rescued is preserved by the merge (rescue is not part of
        // the consensus path today).
        dedup_bam = consensus_out.bam
            .map { meta, bam, bai -> [meta + [bam_dedup: 'xoos_consensus'], bam, bai] }
        dedup_report = consensus_out.report
        dedup_read_collapser_metrics_dir = consensus_out.read_collapser_metrics_dir
        ch_versions = ch_versions.mix(consensus_out.versions)
    }

    emit:
    bam                        = dedup_bam // channel: [ val(meta), path(bam), path(bai) ]
    report                     = dedup_report // channel: [ val(meta), path(json_report) ]
    read_collapser_metrics_dir = dedup_read_collapser_metrics_dir // channel: [ val(meta), path(metrics_dir) ]
    versions                   = ch_versions
}
