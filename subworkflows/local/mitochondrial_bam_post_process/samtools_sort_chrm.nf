include { SAMTOOLS_SORT as SAMTOOLS_SORT_NFCORE } from '../../../modules/nf-core/samtools/sort/main'

workflow SAMTOOLS_SORT_CHRM {
    take:
    ch_with_rg_bam // channel: [ val(meta), path(with_rg_bam) ]

    main:
    SAMTOOLS_SORT_NFCORE(
        ch_with_rg_bam,
        [[], [], []],
        'bai',
    )

    ch_sorted_bam = SAMTOOLS_SORT_NFCORE.out.bam
        .join(SAMTOOLS_SORT_NFCORE.out.index, by: 0)

    emit:
    bam = ch_sorted_bam // channel: [ val(meta), path(chrM_bam), path(chrM_bai) ]
}
