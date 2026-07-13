include { MINIMAP2_ALIGN } from '../../../modules/local/minimap2/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow ALIGN_MINIMAP2 {
    take:
    ch_fastq // [ val(meta), path(reads) ]
    ch_fasta // [ val(meta), path(fasta), path(fai) ]

    main:

    // Prepare reference channel for Minimap2 module
    ch_reference = ch_fasta.map { meta, fasta, _fai -> [meta + [id: 'reference'], fasta] }

    // Align with Minimap2
    mm_out = MINIMAP2_ALIGN(
        ch_fastq.map { meta, reads -> [meta + [single_end: true], reads] },
        ch_reference,
        true
    )

    // Index BAM files
    samtools_index_out = SAMTOOLS_INDEX(mm_out.bam)

    // Combine BAM and BAI into standard output format
    bams = mm_out.bam
        .join(samtools_index_out.index)
        .map { meta, bam, bai -> [meta, bam, bai] }

    emit:
    bams     = bams // channel: [ val(meta), path(bam), path(bai) ]
    versions = Channel.empty() // versions are collected via topic channels
}
