include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_SPLITINTERVALS  } from '../../../modules/nf-core/gatk4/splitintervals/main'
include { GATK4_MERGEVCFS       } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { FAI_TO_BED            } from '../../../modules/local/fai_to_bed/main'

def replaceExtension(path, extension) {
    return file(path).resolveSibling(file(path).baseName + extension)
}

workflow HAPLOTYPECALLER_GATK {
    take:
    ch_bam            // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_reference_dict // channel: [ val(meta), path(dict) ]
    ch_splitintervals_reference      // channel: [ val(meta), path(fasta), path(fai) ] (subset reference for partitioning)
    ch_splitintervals_reference_dict // channel: [ val(meta), path(dict) ] (subset reference dict for partitioning)
    ch_intervals      // channel: [ val(meta), path(bed) ]
    enable_partitioning // val: boolean
    partition_count     // val: integer

    main:

    if (enable_partitioning) {
        // Resolve intervals: use provided BED or fall back to whole-genome BED from FAI.
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
    } else {
        ch_intervals_separated = ch_intervals
            .map { meta, bed -> [(meta instanceof Map ? meta : [:]) + [num_intervals: 1], bed] }
    }

    // Combine BAMs x intervals, adding interval_name to meta
    ch_bam_intervals = ch_bam
        .combine(ch_intervals_separated)
        .map { meta_bam, bam, bai, meta_interval, interval_file ->
            [
                meta_bam + [interval_name: interval_file.baseName, num_intervals: meta_interval.num_intervals],
                bam, bai, interval_file, [],
            ]
        }

    GATK4_HAPLOTYPECALLER(
        ch_bam_intervals,
        ch_reference.map { _meta, fasta, _fai -> [[id: 'fasta'], fasta] },
        ch_reference.map { _meta, _fasta, fai -> [[id: 'fai'], fai] },
        ch_reference_dict.map { _meta, dict -> [[id: 'dict'], dict] },
        channel.value([[], []]),
        channel.value([[], []]),
    )

    // Branch VCFs by whether merging is needed
    hc_vcf = GATK4_HAPLOTYPECALLER.out.vcf
        .map { meta, vcf -> [meta - meta.subMap('interval_name'), vcf] }
        .branch {
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    GATK4_MERGEVCFS(
        hc_vcf.intervals
            .map { meta, vcf -> [groupKey(meta, meta.num_intervals), vcf] }
            .groupTuple(),
        ch_reference_dict.map { _meta, dict -> [[id: 'dict'], dict] },
    )

    // Derive TBI path from VCF for no_intervals branch (TBI output is optional in GATK4_HAPLOTYPECALLER)
    ch_vcfs = Channel.empty()
        .mix(
            GATK4_MERGEVCFS.out.vcf.join(GATK4_MERGEVCFS.out.tbi),
            hc_vcf.no_intervals
                .map { meta, vcf -> [meta, vcf, file("${vcf}.tbi")] },
        )
        .map { meta, vcf, tbi -> [meta - meta.subMap('num_intervals'), vcf, tbi] }

    // Gather BAMs
    hc_bam = GATK4_HAPLOTYPECALLER.out.bam
        .map { meta, bam -> [meta - meta.subMap('interval_name'), bam] }
        .branch {
            intervals:    it[0].num_intervals > 1
            no_intervals: it[0].num_intervals <= 1
        }

    ch_bams = Channel.empty()
        .mix(
            hc_bam.intervals
                .map { meta, bam -> [groupKey(meta, meta.num_intervals), bam] }
                .groupTuple()
                .map { meta, bams -> [meta, bams, bams.collect { replaceExtension(it, '.bai') }] },
            hc_bam.no_intervals
                .map { meta, bam -> [meta, [bam], [replaceExtension(bam, '.bai')]] },
        )
        .map { meta, bams, bais -> [meta - meta.subMap('num_intervals'), bams, bais] }

    emit:
    vcf      = ch_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
    bam      = ch_bams // channel: [ val(meta), [path(bam), ...], [path(bai), ...] ]
}
