/*
 * Scatter-gather subworkflow: partition FASTQ files, align each partition
 * in parallel, then merge the resulting BAMs with sambamba.
 *
 * For giraffe aligners with pangenome personalization, the personalized
 * graph is built once per sample before partitioning, then shared across
 * all partitions.
 *
 * The output channel has the same shape as ALIGN: [ val(meta), path(bam), path(bai) ]
 */

include { ALIGN                    } from './main'
include { PREPARE_GIRAFFE_INPUTS   } from './prepare_giraffe_inputs'
include { VG_GIRAFFE               } from '../../../modules/local/vg/giraffe/main'
include { PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS } from '../../../modules/local/parabricks/giraffe/main'
include { SAMBAMBA_MERGE           } from '../../../modules/local/sambamba/merge/main'
include { partitionFastqs          } from './partition_fastqs'

def isGiraffeAligner(aligner) {
    return aligner in ['vg_giraffe', 'pb_giraffe']
}

workflow PARTITION_ALIGN_MERGE {
    take:
    aligner            // val: aligner name
    partition_count    // val: number of partitions
    ch_fastq           // channel: [ val(meta), path(reads) ]
    ch_reference       // channel: [ [id:'reference'], path(fasta), path(fai) ]
    ch_bwa_idx         // channel: [ [id:'bwa_idx'], path(idx_dir) ]
    ch_giraffe_idx     // channel: [ [id:'giraffe_idx'], path(idx_dir) ]
    ch_ref_paths       // channel: [ path(dict) ] — linearized reference .dict for pbrun giraffe --ref-paths

    main:

    // Partition FASTQs — always emits exactly partition_count items per sample
    ch_all_partitions = partitionFastqs(ch_fastq, partition_count)

    // Separate non-empty partitions (to align) from empty ones (to skip)
    ch_partitioned = ch_all_partitions.filter { _meta, reads -> reads.size() > 0 }
    ch_empty       = ch_all_partitions.filter { _meta, reads -> reads.size() == 0 }

    if (isGiraffeAligner(aligner) && params.enable_pangenome_personalization.toBoolean()) {
        // Giraffe + personalization: build graph once per sample, then partition

        // 1. Build personalized graphs from unpartitioned FASTQs (one per sample)
        giraffe_inputs = PREPARE_GIRAFFE_INPUTS(
            ch_fastq,
            ch_giraffe_idx,
            true,
        ).align_inputs

        // 2. Extract graph components keyed by sample id (for later recombination)
        ch_graphs = giraffe_inputs
            .map { meta, _reads, graph, zipcodes, minimizers, dist ->
                [meta.id, graph, zipcodes, minimizers, dist]
            }

        // 3. Combine each partition with the shared graph using original_id as key
        ch_partitioned_with_graphs = ch_partitioned
            .map { meta, reads -> [meta.original_id, meta, reads] }
            .combine(ch_graphs, by: 0)
            .map { _key, meta, reads, graph, zipcodes, minimizers, dist ->
                [meta, reads, graph, zipcodes, minimizers, dist]
            }

        // 4. Align each partition with the shared personalized graph
        if (aligner == 'vg_giraffe') {
            align_bams = VG_GIRAFFE(ch_partitioned_with_graphs).bam
        } else {
            align_bams = PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS(ch_partitioned_with_graphs, ch_ref_paths).bam
        }
    } else {
        // Non-giraffe or giraffe without personalization: delegate to ALIGN
        align_bams = ALIGN(
            aligner,
            ch_partitioned,
            ch_reference,
            ch_bwa_idx,
            ch_giraffe_idx,
            ch_ref_paths,
        ).bams
    }

    // Gather: group partition BAMs back by sample and merge
    //   - Restore original_id back to id and remove partition keys
    //   - Use groupTuple(size: partition_count) for eager emission since
    //     partitionFastqs always emits exactly partition_count items per sample
    //   - Empty partitions contribute a null BAM placeholder to reach the expected count

    // Create null BAM placeholders for empty partitions so groupTuple gets
    // exactly partition_count items per sample
    ch_empty_bams = ch_empty.map { meta, _reads ->
        [meta, null, null]
    }

    ch_grouped_bams = align_bams
        .mix(ch_empty_bams)
        .map { meta, bam, bai ->
            def clean_meta = meta.clone()
            clean_meta.remove('partition')
            clean_meta.remove('single_end') // Added by ALIGN subworkflows; must be removed so all partitions share identical meta for groupTuple
            clean_meta.id = clean_meta.remove('original_id')
            [clean_meta, bam, bai]
        }
        .groupTuple(size: partition_count as int)
        .map { meta, bams, bais ->
            // Filter out null placeholders from empty partitions
            def real_bams = bams.flatten().findAll { it != null }
            def real_bais = bais.flatten().findAll { it != null }
            [meta, real_bams, real_bais]
        }

    // sambamba merge requires >= 2 inputs; pass single-BAM samples through
    ch_to_merge  = ch_grouped_bams.filter { _meta, bams, _bais -> bams.size() > 1 }
    ch_single    = ch_grouped_bams.filter { _meta, bams, _bais -> bams.size() == 1 }
        .map { meta, bams, bais -> [meta, bams[0], bais[0]] }

    bams = SAMBAMBA_MERGE(ch_to_merge).bam.mix(ch_single)

    emit:
    bams // channel: [ val(meta), path(bam), path(bai) ]
}
