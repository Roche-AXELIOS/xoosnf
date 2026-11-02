include { MARKDUP               } from './markdup'
include { CONSENSUS_AND_REALIGN } from './consensus'

workflow DEDUP {
    take:
    dedup_strategy // val: 'markdup' or 'consensus'
    aligner // val: 'bwa', 'minimap2', 'vg_giraffe', 'pb_fq2bam', 'pb_minimap2', or 'pb_giraffe'
    enable_secondary_rescue // val: true if this run needs secondary-alignment rescue
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
    dedup_giraffe_graphs = channel.empty()
    dedup_report = channel.empty()
    dedup_read_collapser_metrics_dir = channel.empty()
    ch_markdup_input_bams = channel.empty()
    ch_consensus_input_bams = channel.empty()
    ch_externally_deduped_bams = channel.empty()
    if (dedup_strategy == 'markdup') {
        // Skip marking only when the input already satisfies this run's
        // requested markdup/rescue state. sas_markdup is deliberately not
        // trusted, so it always continues to XOOS read_collapser MARKDUP.
        ch_externally_deduped_bams = ch_bam.filter { meta, _bam, _bai ->
            meta.bam_dedup_type == 'xoos_markdup' &&
                (!enable_secondary_rescue || meta.bam_secondary_rescued)
        }
        ch_markdup_input_bams = ch_bam.filter { meta, _bam, _bai ->
            !(meta.bam_dedup_type == 'xoos_markdup' &&
                (!enable_secondary_rescue || meta.bam_secondary_rescued))
        }

        markdup_out = MARKDUP(ch_markdup_input_bams, enable_secondary_rescue)
        dedup_bam = markdup_out.markdup_bam
            .map { meta, bam, bai -> [meta + [bam_dedup_type: 'xoos_markdup', bam_secondary_rescued: meta.bam_secondary_rescued ?: false], bam, bai] }
            .mix(ch_externally_deduped_bams)
        dedup_report = markdup_out.report
        dedup_read_collapser_metrics_dir = markdup_out.read_collapser_metrics_dir
        ch_versions = ch_versions.mix(markdup_out.versions)
    }
    else if (dedup_strategy == 'consensus') {
        // Secondary-alignment rescue is only wired into the markdup path. If the
        // run configuration requests rescue while using the consensus strategy,
        // it will not be applied — warn so this silent no-op is visible.
        if (enable_secondary_rescue) {
            log.warn("enable_secondary_rescue is enabled but dedup_strategy is 'consensus'; secondary-alignment rescue is only applied in the markdup path and will be skipped.")
        }
        ch_externally_deduped_bams = ch_bam.filter { meta, _bam, _bai ->
            meta.bam_dedup_type == 'xoos_consensus'
        }
        ch_consensus_input_bams = ch_bam.filter { meta, _bam, _bai ->
            meta.bam_dedup_type != 'xoos_consensus'
        }

        consensus_out = CONSENSUS_AND_REALIGN(
            ch_consensus_input_bams,
            ch_bed,
            aligner,
            ch_fasta,
            ch_bwa_idx,
            ch_giraffe_idx,
            ch_ref_paths,
        )
        dedup_bam = consensus_out.bam
            .map { meta, bam, bai -> [meta + [bam_dedup_type: 'xoos_consensus'], bam, bai] }
            .mix(ch_externally_deduped_bams)
        dedup_giraffe_graphs = consensus_out.giraffe_graphs
        dedup_report = consensus_out.report
        dedup_read_collapser_metrics_dir = consensus_out.read_collapser_metrics_dir
        ch_versions = ch_versions.mix(consensus_out.versions)
    }

    emit:
    bam                        = dedup_bam // channel: [ val(meta), path(bam), path(bai) ]
    giraffe_graphs             = dedup_giraffe_graphs // channel: [ val(meta), path(gbz), path(zipcodes), path(min), path(dist) ]
    report                     = dedup_report // channel: [ val(meta), path(json_report) ]
    read_collapser_metrics_dir = dedup_read_collapser_metrics_dir // channel: [ val(meta), path(metrics_dir) ]
    versions                   = ch_versions
}
