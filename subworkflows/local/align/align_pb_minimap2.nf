include { PARABRICKS_MINIMAP2 } from '../../../modules/local/parabricks/minimap2/main'

workflow ALIGN_PB_MINIMAP2 {
    take:
    ch_fastq // channel: [ val(meta), path(reads) ]
    ch_fasta // channel: [ val(meta), path(fasta), path(fai) ]

    main:
    // Prepare reference channel for Parabricks module
    ch_reference = ch_fasta.map { meta, fasta, _fai -> [meta + [id: 'reference'], fasta] }
    empty_paths = channel.value([[], []])

    // Convert XOOS meta format to nf-core format
    nfcore_fastqs = ch_fastq.map { meta, reads -> [meta + [single_end: true], reads] }

    // Align with Parabricks minimap2 (GPU-accelerated)
    pb_out = PARABRICKS_MINIMAP2(
        nfcore_fastqs,
        ch_reference,
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
