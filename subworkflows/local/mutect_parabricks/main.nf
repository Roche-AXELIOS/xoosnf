include { GATK4_SPLITINTERVALS    } from '../../../modules/nf-core/gatk4/splitintervals/main'
include { GATK4_MERGEVCFS         } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEMUTECTSTATS  } from '../../../modules/local/gatk4/mergemutectstats/main'
include { FAI_TO_BED              } from '../../../modules/local/fai_to_bed/main'
include { PARABRICKS_MUTECTCALLER } from '../../../modules/local/parabricks/mutectcaller/main'

workflow MUTECT_PARABRICKS {
    take:
    ch_paired_bams    // channel: [ val(meta), [path(bam), ...], [path(bai), ..] ]
    ch_reference      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_reference_dict // channel: [ val(meta), path(dict) ]
    ch_splitintervals_reference      // channel: [ val(meta), path(fasta), path(fai) ] (subset reference for partitioning)
    ch_splitintervals_reference_dict // channel: [ val(meta), path(dict) ] (subset reference dict for partitioning)
    ch_intervals      // channel: [ val(meta), path(intervals) ]
    enable_mutect_stats // val: boolean
    enable_partitioning // val: boolean
    partition_count     // val: integer

    main:
    // Collect versions (also published via topic channels)
    ch_versions = Channel.empty()

    if (enable_partitioning) {
        // Resolve intervals: use provided BED or fall back to whole-genome BED from FAI.
        // The subset reference drives partitioning so that variant calling is split only
        // across the subset's chromosomes; the callers below still use the full reference.
        ch_whole_genome_bed = FAI_TO_BED(ch_splitintervals_reference.map { meta, _fasta, fai -> [meta, fai] })
        ch_resolved_intervals = ch_intervals
            .combine(ch_whole_genome_bed.bed)
            .map { _meta, bed, meta2, whole_genome_bed ->
                [meta2 + [variant_calling_partitions: partition_count], bed ?: whole_genome_bed]
            }

        GATK4_SPLITINTERVALS(
            ch_resolved_intervals,
            ch_splitintervals_reference.map { _meta, fasta, _fai -> [[id: 'fasta'], fasta] },
            ch_splitintervals_reference.map { _meta, _fasta, fai -> [[id: 'fai'], fai] },
            ch_splitintervals_reference_dict.map { _meta, dict -> [[id: 'dict'], dict] },
        )

        // Derive num_intervals from actual SPLITINTERVALS output count (may differ from partition_count)
        ch_intervals_separated = GATK4_SPLITINTERVALS.out.split_intervals.flatMap { meta, intervals ->
            def interval_list = intervals instanceof Collection ? intervals.toList() : [intervals]
            interval_list.collect { interval_file -> [meta + [num_intervals: interval_list.size()], interval_file] }
        }
    } else {
        ch_intervals_separated = ch_intervals
            .map { meta, bed -> [(meta instanceof Map ? meta : [:]) + [num_intervals: 1], bed] }
    }

    // Expand paired BAMs into individual tumor/normal components
    ch_expanded_bams = ch_paired_bams
        .map { meta, bams, bais ->
            bams.size() == 2
                ? [meta, bams[0], bais[0], bams[1], bais[1]]
                : [meta, bams[0], bais[0], [], []]
        }

    // Combine BAMs x intervals, adding interval_name to meta
    ch_bam_intervals = ch_expanded_bams
        .combine(ch_intervals_separated)
        .map { meta, bam1, bai1, bam2, bai2, meta_interval, interval_file ->
            [
                meta + [interval_name: interval_file.baseName, num_intervals: meta_interval.num_intervals],
                bam1, bai1, bam2, bai2, interval_file,
            ]
        }

    PARABRICKS_MUTECTCALLER(
        ch_bam_intervals,
        ch_reference.map { meta, fasta, _fai -> [meta, fasta] },
        [],
        [],
    )

    // Branch VCFs by whether merging is needed (Parabricks emits [meta, vcf, tbi])
    mutect_vcf = PARABRICKS_MUTECTCALLER.out.vcf
        .map { meta, vcf, _tbi -> [meta - meta.subMap('interval_name'), vcf] }
        .branch {
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    GATK4_MERGEVCFS(
        mutect_vcf.intervals
            .map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }
            .groupTuple(),
        ch_reference_dict.map { _meta, dict -> [[id: 'dict'], dict] },
    )

    ch_vcfs = Channel.empty()
        .mix(
            GATK4_MERGEVCFS.out.vcf.join(GATK4_MERGEVCFS.out.tbi),
            PARABRICKS_MUTECTCALLER.out.vcf
                .map { meta, vcf, tbi -> [meta - meta.subMap('interval_name'), vcf, tbi] }
                .filter { meta, _vcf, _tbi -> meta.num_intervals <= 1 },
        )
        .map { meta, vcf, tbi -> [meta - meta.subMap('num_intervals'), vcf, tbi] }

    // Gather BAMs (Parabricks emits [meta, bam, bai])
    mutect_bam = PARABRICKS_MUTECTCALLER.out.bam
        .map { meta, bam, bai -> [meta - meta.subMap('interval_name'), bam, bai] }
        .branch {
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    ch_bams_out = Channel.empty()
        .mix(
            mutect_bam.intervals
                .map { meta, bam, bai -> [groupKey(meta, meta.num_intervals), bam, bai] }
                .groupTuple()
                .map { meta, bams, bais -> [meta, bams, bais] },
            mutect_bam.no_intervals
                .map { meta, bam, bai -> [meta, [bam], [bai]] },
        )
        .map { meta, bams, bais -> [meta - meta.subMap('num_intervals'), bams, bais] }

    // Gather stats (optional)
    ch_merged_stats = channel.empty()
    if (enable_mutect_stats) {
        mutect_stats = PARABRICKS_MUTECTCALLER.out.stats
            .map { meta, stats -> [meta - meta.subMap('interval_name'), stats] }
            .branch {
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }

        GATK4_MERGEMUTECTSTATS(
            mutect_stats.intervals
                .map { meta, stats -> [groupKey(meta, meta.num_intervals), stats] }
                .groupTuple()
                .map { meta, stats -> [meta, stats.flatten()] }
        )

        ch_merged_stats = Channel.empty()
            .mix(GATK4_MERGEMUTECTSTATS.out.stats, mutect_stats.no_intervals)
            .map { meta, stats -> [meta - meta.subMap('num_intervals'), stats] }
    }

    emit:
    vcf      = ch_vcfs         // channel: [ val(meta), path(vcf), path(tbi) ]
    bam      = ch_bams_out     // channel: [ val(meta), [path(bam), ...], [path(bai), ...] ]
    stats    = ch_merged_stats // channel: [ val(meta), path(stats) ]
    versions = ch_versions
}
