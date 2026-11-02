/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DOWNSAMPLED BWA METRICS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downsamples raw consensus FASTQs, re-aligns them with the configured
    aligner (bwa or pb_fq2bam, selected via the `aligner` workflow input),
    and runs the full alignment metrics suite on the re-aligned BAM. Used for
    cross-validation with the R&D pipeline (parity check).

    This runs in parallel with the primary alignment path and is gated by
    enable_downsampled_bwa_metrics (off by default).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SUBSAMPLE_FASTQ_RANDOM } from '../../../modules/local/subsample/fastq/seqkit_random/main'
include { ALIGN_BWA              } from '../align/align_bwa'
include { ALIGN_PB_FQ2BAM        } from '../align/align_pb_fq2bam'
include { METRICS as BWA_METRICS } from '../metrics/main'

workflow DOWNSAMPLED_BWA_METRICS {
    take:
    ch_fastq                            // channel: [ val(meta), path(reads) ]
    ch_reference                        // channel: [ val(meta), path(fasta), path(fai) ]
    ch_bwa_idx                          // channel: [ val(meta), path(index) ]
    ch_autosomes_bed                    // channel: [ val(meta), path(bed) ]
    ch_target_bed                       // channel: [ val(meta), path(bed) ]
    ch_capture_bed                      // channel: [ val(meta), path(bed) ]
    ch_high_confidence_bed_dir          // channel: [ val(meta), path(dir) ]
    ch_fallback_high_confidence_bed     // channel: [ val(meta), path(bed) ]
    ch_reference_for_metrics                 // channel: [ val(meta), path(fasta), path(fai) ]
    enable_alignment_metrics_autosomes  // boolean
    enable_alignment_metrics_targets    // boolean
    enable_alignment_metrics_high_confidence // boolean
    enable_gc_bias_metrics              // boolean
    enable_hs_metrics                   // boolean
    downsampled_bwa_metrics_read_count  // integer: number of reads to downsample to
    downsampled_bwa_metrics_aligner     // string: aligner to use ('bwa' or 'pb_fq2bam')
    additional_bed_files                // list: [ [region_name, bed] ]

    main:
    ch_versions = Channel.empty()

    // Step 1: Downsample FASTQs with seqkit sample
    // Input: [ meta, read_count, fastqs ]
    ch_fastq_with_read_count = ch_fastq.map { meta, reads ->
        [meta, downsampled_bwa_metrics_read_count, reads]
    }
    SUBSAMPLE_FASTQ_RANDOM(ch_fastq_with_read_count)
    ch_versions = ch_versions.mix(SUBSAMPLE_FASTQ_RANDOM.out.versions)

    // Step 2: Align downsampled FASTQs
    // Select aligner based on the aligner input
    // Both ALIGN_BWA and ALIGN_PB_FQ2BAM share the same interface:
    //   take: ch_fastq, ch_fasta, ch_bwa_idx
    //   emit: bams [ val(meta), path(bam), path(bai) ]
    ch_aligned_bams = channel.empty()
    if (downsampled_bwa_metrics_aligner == 'bwa') {
        ALIGN_BWA(
            SUBSAMPLE_FASTQ_RANDOM.out.fastq.map { meta, fq -> [meta, [fq]] },
            ch_reference,
            ch_bwa_idx,
        )
        ch_aligned_bams = ALIGN_BWA.out.bams
        ch_versions = ch_versions.mix(ALIGN_BWA.out.versions)
    } else if (downsampled_bwa_metrics_aligner == 'pb_fq2bam') {
        // Default: pb_fq2bam (GPU-accelerated)
        ALIGN_PB_FQ2BAM(
            SUBSAMPLE_FASTQ_RANDOM.out.fastq.map { meta, fq -> [meta, [fq]] },
            ch_reference,
            ch_bwa_idx,
        )
        ch_aligned_bams = ALIGN_PB_FQ2BAM.out.bams
        ch_versions = ch_versions.mix(ALIGN_PB_FQ2BAM.out.versions)
    } else {
        error "Invalid value for aligner: ${downsampled_bwa_metrics_aligner}. Must be one of: bwa, pb_fq2bam"
    }

    // Step 3: Run full metrics suite on aligned BAM
    BWA_METRICS(
        'downsampled_metrics',
        ch_aligned_bams,
        ch_autosomes_bed,
        ch_target_bed,
        ch_capture_bed,
        ch_high_confidence_bed_dir,
        ch_fallback_high_confidence_bed,
        ch_reference_for_metrics,
        enable_alignment_metrics_autosomes,
        enable_alignment_metrics_targets,
        enable_alignment_metrics_high_confidence,
        enable_gc_bias_metrics,
        enable_hs_metrics,
        downsampled_bwa_metrics_aligner,
        additional_bed_files,
    )
    ch_versions = ch_versions.mix(BWA_METRICS.out.versions)

    emit:
    report                  = BWA_METRICS.out.report                  // channel: [ val(meta), path(report) ]
    alignment_metrics_files = BWA_METRICS.out.alignment_metrics_files // channel: [ val(meta), path(read_metrics), path(coverage), path(accuracy) ]
    gc_bias_metrics_files   = BWA_METRICS.out.gc_bias_metrics_files   // channel: [ val(meta), path(metrics) ]
    hs_metrics_files        = BWA_METRICS.out.hs_metrics_files        // channel: [ val(meta), path(metrics) ]
    custom_multiqc_files    = BWA_METRICS.out.custom_multiqc_files
    multiqc_files           = BWA_METRICS.out.multiqc_files
    versions                = ch_versions
}
