include { MUTECT_GATK as MUTECT_GATK_TN                              } from '../mutect_gatk'
include { MUTECT_PARABRICKS as MUTECT_PARABRICKS_TN                  } from '../mutect_parabricks'
include { XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_TUMOR_NORMAL_WGS } from '../../../modules/local/xoos_small_variant_caller/filter_variants_tumor_normal_wgs/main'

workflow CALL_SOMATIC_SMALL_VARIANTS_TN {
    take:
    ch_paired_bams // channel: [meta, bam_tumor, bai_tumor, bam_normal, bai_normal]
    ch_reference // channel: [meta, fasta, fai]
    ch_reference_dict // channel: [meta, dict]
    ch_splitintervals_reference // channel: [meta, fasta, fai] (subset reference for partitioning)
    ch_splitintervals_reference_dict // channel: [meta, dict] (subset reference dict for partitioning)
    ch_pop_af_vcf_with_tbi // channel: [ val(meta), path(vcf.gz), path(vcf.gz.tbi) ]
    ch_svc_model // channel: [meta, model]
    ch_intervals_bed // channel: [meta, bed]
    use_gatk_mutect2 // val
    enable_mutect_stats // val: boolean
    enable_partitioning // val: boolean
    partition_count // val: integer
    aligner // val: string

    main:
    ch_mutect_vcf = channel.empty()
    ch_mutect_bam = channel.empty()
    if (use_gatk_mutect2) {
        MUTECT_GATK_TN(
            ch_paired_bams,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_intervals_bed,
            enable_mutect_stats,
            enable_partitioning,
            partition_count,
        )
        ch_mutect_vcf = MUTECT_GATK_TN.out.vcf
        ch_mutect_bam = MUTECT_GATK_TN.out.bam
    }
    else {
        MUTECT_PARABRICKS_TN(
            ch_paired_bams,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_intervals_bed,
            enable_mutect_stats,
            enable_partitioning,
            partition_count,
        )
        ch_mutect_vcf = MUTECT_PARABRICKS_TN.out.vcf
        ch_mutect_bam = MUTECT_PARABRICKS_TN.out.bam
    }
    XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_TUMOR_NORMAL_WGS(
        ch_mutect_bam.join(ch_mutect_vcf),
        ch_reference,
        ch_svc_model,
        ch_pop_af_vcf_with_tbi,
        ch_intervals_bed,
        aligner,
    )

    emit:
    raw_vcf     = ch_mutect_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    vcf         = XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_TUMOR_NORMAL_WGS.out.vcf
}
