include { PARABRICKS_HAPLOTYPECALLER } from '../../../modules/local/parabricks/haplotypecaller/main'
include { HTSLIB_BGZIPTABIX           } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { GATK4_SPLITINTERVALS       } from '../../../modules/nf-core/gatk4/splitintervals/main'
include { GATK4_MERGEVCFS            } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { FAI_TO_BED                 } from '../../../modules/local/fai_to_bed/main'

workflow HAPLOTYPECALLER_PARABRICKS {
    take:
    bams              // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_reference_dict // channel: [ val(meta), path(dict) ]
    ch_splitintervals_reference      // channel: [ val(meta), path(fasta), path(fai) ] (subset reference for partitioning)
    ch_splitintervals_reference_dict // channel: [ val(meta), path(dict) ] (subset reference dict for partitioning)
    ch_intervals      // channel: [ val(meta), path(intervals) ]
    enable_partitioning // val: boolean
    partition_count     // val: integer

    main:

    if (enable_partitioning) {
        // Resolve intervals: use provided BED or fall back to whole-genome BED from FAI
        // (GATK4_SPLITINTERVALS requires an interval file).
        // The subset reference drives partitioning so that variant calling is split only
        // across the subset's chromosomes; the callers below still use the full reference.
        ch_whole_genome_bed = FAI_TO_BED(ch_splitintervals_reference.map { meta, _fasta, fai -> [meta, fai] })
        ch_resolved_intervals = ch_intervals.combine(ch_whole_genome_bed.bed)
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

        // Combine BAMs × intervals, adding interval_name to meta
        ch_bam_intervals = bams
            .combine(ch_intervals_separated)
            .map { meta_bam, bam, bai, meta_interval, interval_file ->
                [
                    meta_bam + [interval_name: interval_file.baseName, num_intervals: meta_interval.num_intervals],
                    bam, bai, interval_file,
                ]
            }
    } else {
        // No partitioning — pass original intervals directly (may be empty)
        ch_bam_intervals = bams
            .combine(ch_intervals)
            .map { meta_bam, bam, bai, _meta_interval, interval_file ->
                [
                    meta_bam + [num_intervals: 1],
                    bam, bai, interval_file,
                ]
            }
    }

    PARABRICKS_HAPLOTYPECALLER(
        ch_bam_intervals,
        ch_reference,
    )

    // Parabricks outputs uncompressed VCF — bgzip+index
    HTSLIB_BGZIPTABIX(
        PARABRICKS_HAPLOTYPECALLER.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
        "compress",
        true,
        "vcf",
    )

    // Branch VCFs by whether merging is needed
    pb_vcf = HTSLIB_BGZIPTABIX.out.output
        .map { meta, vcf -> [meta - meta.subMap('interval_name'), vcf] }
        .branch {
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    GATK4_MERGEVCFS(
        pb_vcf.intervals
            .map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }
            .groupTuple(),
        ch_reference_dict.map { _meta, dict -> [[id: 'dict'], dict] },
    )

    ch_vcfs = Channel.empty()
        .mix(
            GATK4_MERGEVCFS.out.vcf.join(GATK4_MERGEVCFS.out.tbi),
            HTSLIB_BGZIPTABIX.out.output.join(HTSLIB_BGZIPTABIX.out.index)
                .map { meta, vcf, tbi -> [meta - meta.subMap('interval_name'), vcf, tbi] }
                .filter { meta, _vcf, _tbi -> meta.num_intervals <= 1 },
        )
        .map { meta, vcf, tbi -> [meta - meta.subMap('num_intervals'), vcf, tbi] }

    // Gather BAMs (Parabricks emits [meta, bam, bai])
    pb_bam = PARABRICKS_HAPLOTYPECALLER.out.bam
        .map { meta, bam, bai -> [meta - meta.subMap('interval_name'), bam, bai] }
        .branch {
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    // Flatten after groupTuple: multi-GPU Parabricks emits multiple BAMs per
    // partition (one per GPU), so groupTuple produces nested lists that must
    // be flattened into a single flat list for downstream path() inputs.
    ch_bams_out = Channel.empty()
        .mix(
            pb_bam.intervals
                .map { meta, bam, bai -> [groupKey(meta, meta.num_intervals), bam, bai] }
                .groupTuple()
                .map { meta, merged_bams, merged_bais -> [meta, merged_bams.flatten(), merged_bais.flatten()] },
            pb_bam.no_intervals
                .map { meta, bam, bai ->
                    [meta, bam instanceof List ? bam : [bam], bai instanceof List ? bai : [bai]]
                },
        )
        .map { meta, merged_bams, merged_bais -> [meta - meta.subMap('num_intervals'), merged_bams, merged_bais] }

    emit:
    vcf      = ch_vcfs     // channel: [ val(meta), path(vcf), path(tbi) ]
    bam      = ch_bams_out // channel: [ val(meta), [path(bam), ...], [path(bai), ...] ]
}
