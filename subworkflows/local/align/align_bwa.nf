include { BWA_MEM_SINGLE_END_MULTIPLE_FASTQS } from '../../../modules/local/bwa/mem/main'
include { SAMTOOLS_INDEX                     } from '../../../modules/nf-core/samtools/index/main'

workflow ALIGN_BWA {
    take:
    ch_fastq // [ val(meta), path(reads) ]
    ch_fasta // [ val(meta), path(fasta), path(fai) ]
    ch_bwa_idx // [ val(meta), path(index) ]

    main:

    // Prepare reference channels for BWA module (nf-core BWA_MEM expects specific format)
    ch_reference = ch_fasta.map { meta, fasta, _fai -> [meta + [id: 'reference'], fasta] }
    ch_ref_index = ch_bwa_idx.map { meta, index -> [meta + [id: 'reference'], index] }

    // Align with BWA-MEM
    bwa_out = BWA_MEM_SINGLE_END_MULTIPLE_FASTQS(ch_fastq.map { meta, reads -> [meta + [single_end: true], reads] }, ch_ref_index, ch_reference, true)

    // Index BAM files
    samtools_index_out = SAMTOOLS_INDEX(bwa_out.bam)

    // Combine BAM and BAI into standard output format
    bams = bwa_out.bam
        .join(samtools_index_out.index)
        .map { meta, bam, bai -> [meta, bam, bai] }

    emit:
    bams     = bams // channel: [ val(meta), path(bam), path(bai) ]
    versions = Channel.empty() // versions are collected via topic channels
}
