include { PARABRICKS_FQ2BAM_SINGLE_END_MULTIPLE_FASTQS } from '../../../modules/local/parabricks/fq2bam/main'

workflow ALIGN_PB_FQ2BAM {
    take:
    ch_fastq // channel: [ val(meta), path(reads) ]
    ch_fasta // channel: [ val(meta), path(fasta), path(fai) ]
    ch_bwa_idx // channel: [ val(meta), path(index) ]

    main:
    // Prepare reference channels for Parabricks module
    ch_reference = ch_fasta.map { meta, fasta, _fai -> [meta + [id: 'reference'], fasta] }
    ref_index_ch = ch_bwa_idx.map { meta, index -> [meta + [id: 'reference'], index] }
    empty_paths = channel.value([[], []])

    // Convert XOOS meta format to nf-core format
    nfcore_fastqs = ch_fastq.map { meta, reads -> [meta + [single_end: true], reads] }

    // Align with Parabricks fq2bam (GPU-accelerated)
    pb_out = PARABRICKS_FQ2BAM_SINGLE_END_MULTIPLE_FASTQS(
        nfcore_fastqs,
        ch_reference,
        ref_index_ch,
        empty_paths,
        empty_paths,
    )

    // Combine BAM and BAI
    bams = pb_out.bam
        .join(pb_out.bai)
        .map { meta, bam, bai -> [meta, bam, bai] }

    emit:
    bams     = bams // channel: [ val(meta), path(bam), path(bai) ]
    versions = Channel.empty() // versions are collected via topic channels
}
