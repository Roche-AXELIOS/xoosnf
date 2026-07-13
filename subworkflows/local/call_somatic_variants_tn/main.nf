include { CALL_SOMATIC_SMALL_VARIANTS_TN                                     } from './call_somatic_small_variants_tn'
include { EVALUATE_SMALL_VARIANT_CALLS as EVALUATE_SOMATIC_SMALL_VARIANTS_TN } from '../evaluate_small_variant_calls/main'

include { BCFTOOLS_CONCAT                                                    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_RENAME_TUMOR_SAMPLE         } from '../../../modules/nf-core/bcftools/reheader/main'
include { XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS                           } from '../../../modules/local/xoos_copy_number_caller/tumor_normal_wgs/main'
include { BCFTOOLS_INDEX                                                     } from '../../../modules/nf-core/bcftools/index/main'

workflow CALL_SOMATIC_VARIANTS_TN {
    take:
    ch_tumor_normal_bams
    ch_reference
    ch_reference_dict
    ch_splitintervals_reference
    ch_splitintervals_reference_dict
    ch_intervals
    ch_pop_af_vcf_with_tbi
    ch_svc_model_tumor_normal
    ch_truth_vcf_dir_tumor_normal
    ch_eval_bed_dir_tumor_normal
    ch_mappability_bw
    ch_cnv_seed_file
    ch_cnv_blocklist_bed
    enable_small_variant_caller
    use_gatk_mutect2
    enable_mutect_stats
    enable_vcfeval
    enable_copy_number_caller
    enable_merge_vcf
    reference_sdf
    enable_partitioning
    partition_count
    aligner

    main:

    ch_all_vcfs = channel.empty()
    ch_raw_vc_vcf = channel.empty()
    ch_vcfeval_report = channel.empty()
    ch_vcfeval_output = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_copy_number_caller_output = channel.empty()
    ch_tumor_normal_filtered_vcf = channel.empty()

    // Step 1: Call small variants (Mutect2 + pretrained filtering)
    // Note: Somatic TN copy number calling requires the small variant VCF as input
    if (enable_small_variant_caller || enable_copy_number_caller) {
        CALL_SOMATIC_SMALL_VARIANTS_TN(
            ch_tumor_normal_bams,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_pop_af_vcf_with_tbi,
            ch_svc_model_tumor_normal,
            ch_intervals,
            use_gatk_mutect2,
            enable_mutect_stats,
            enable_partitioning,
            partition_count,
            aligner,
        )
        ch_raw_vc_vcf = CALL_SOMATIC_SMALL_VARIANTS_TN.out.raw_vcf
        ch_tumor_normal_filtered_vcf = CALL_SOMATIC_SMALL_VARIANTS_TN.out.vcf

        // Add pretrained VCFs to merged output
        ch_all_vcfs = ch_all_vcfs.mix(ch_tumor_normal_filtered_vcf)

        // Step 2: Evaluate somatic small variant calls with RTG vcfeval, if enabled
        if (enable_vcfeval) {
            if (enable_small_variant_caller) {
                // Filter out normal sample from the VCF input to vcfeval
                // and rename the tumor sample to match the truth VCF sample name
                // This is because vcfeval expects a single sample column with the same name
                // as the tumor-normal group name in both the query and truth VCFs
                ch_sample_names = ch_tumor_normal_filtered_vcf.collectFile { meta, _vcf, _tbi ->
                    ["${meta.id}.txt", "${meta.tumor_id} ${meta.id}" + '\n']
                }
                ch_bcftools_reheader_input = ch_tumor_normal_filtered_vcf
                    .combine(ch_sample_names)
                    .filter { meta, _vcf, _tbi, sample_file -> sample_file.baseName == meta.id }
                    .map { meta, vcf, _tbi, sample_file -> [meta, vcf, [], sample_file] }
                BCFTOOLS_REHEADER_RENAME_TUMOR_SAMPLE(
                    ch_bcftools_reheader_input,
                    channel.value([[], []]),
                )
                ch_vcfeval_input = BCFTOOLS_REHEADER_RENAME_TUMOR_SAMPLE.out.vcf.join(BCFTOOLS_REHEADER_RENAME_TUMOR_SAMPLE.out.index)

                ch_vcfeval_report = EVALUATE_SOMATIC_SMALL_VARIANTS_TN(
                    ch_vcfeval_input.map { meta, vcf, tbi -> [meta + [svc_model_type: 'pretrained'], vcf, tbi] },
                    ch_reference,
                    ch_truth_vcf_dir_tumor_normal,
                    ch_eval_bed_dir_tumor_normal,
                    reference_sdf,
                ).report.concat(ch_vcfeval_report)
                ch_multiqc_files = EVALUATE_SOMATIC_SMALL_VARIANTS_TN.out.custom_multiqc_files.concat(ch_multiqc_files)
                ch_vcfeval_output = EVALUATE_SOMATIC_SMALL_VARIANTS_TN.out.vcfeval_output.concat(ch_vcfeval_output)
            }
        }

        // Step 3: Call copy number variants with XOOS Copy Number Caller, if enabled
        if (enable_copy_number_caller) {
            ch_cnc_input = ch_tumor_normal_bams
                .map { meta, tumor_normal_bams, tumor_normal_bam_indices -> [meta, tumor_normal_bams[0], tumor_normal_bam_indices[0], tumor_normal_bams[1], tumor_normal_bam_indices[1]] }
                .join(ch_tumor_normal_filtered_vcf)
            ch_copy_number_vcf = XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS(
                ch_cnc_input,
                ch_reference,
                ch_mappability_bw,
                ch_cnv_seed_file,
                ch_cnv_blocklist_bed,
            ).vcf
            ch_copy_number_caller_output = XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS.out.results_dir
            ch_copy_number_vcf_index = BCFTOOLS_INDEX(ch_copy_number_vcf).index
            ch_all_vcfs = ch_copy_number_vcf
                .join(ch_copy_number_vcf_index)
                .mix(ch_all_vcfs)
        }
    }

    ch_merged_vcf = channel.empty()
    if (enable_merge_vcf) {
        // Count enabled callers so groupTuple knows the expected group size per sample.
        // SVC is always run when we reach this point (guard: enable_small_variant_caller || enable_copy_number_caller).
        def num_callers = 1 + (enable_copy_number_caller ? 1 : 0)
        ch_merged_vcf = BCFTOOLS_CONCAT(
            ch_all_vcfs
                .map { meta, vcf, tbi -> [groupKey(meta, num_callers), vcf, tbi] }
                .groupTuple()
        ).vcf
        ch_merged_vcf = ch_merged_vcf.join(BCFTOOLS_CONCAT.out.index)
    } else {
        // Emit individual VCFs without merging
        ch_merged_vcf = ch_all_vcfs.map { meta, vcf, tbi -> [meta, vcf, tbi] }
    }

    emit:
    // When enable_merge_vcf is true: one [meta, vcf, tbi] tuple per sample (merged).
    // When enable_merge_vcf is false: one [meta, vcf, tbi] tuple per caller per sample (unmerged).
    merged_vcf                   = ch_merged_vcf
    raw_vc_vcf                   = ch_raw_vc_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    svc_vcf                      = ch_tumor_normal_filtered_vcf // [ meta, path(vcf), path(tbi) ]
    small_variant_vcfeval_report = ch_vcfeval_report
    small_variant_vcfeval_output = ch_vcfeval_output // channel: [ val(meta), path(vcfeval_dir) ]
    copy_number_caller_output    = ch_copy_number_caller_output
    custom_multiqc_files         = ch_multiqc_files
}
