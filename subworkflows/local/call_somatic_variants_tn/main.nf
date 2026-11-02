include { CALL_SOMATIC_SMALL_VARIANTS_TN                                     } from './call_somatic_small_variants_tn'
include { EVALUATE_SMALL_VARIANT_CALLS as EVALUATE_SOMATIC_SMALL_VARIANTS_TN } from '../evaluate_small_variant_calls/main'
include { resolveTruthSampleName                                              } from '../utils_samplesheet/main'
include { RESOLVE_TUMOR_NORMAL_BAM_SAMPLE_NAMES                               } from './resolve_tumor_normal_bam_sample_names'
include { REPLACE_VCF_IN_RESULTS as REPLACE_CNC_TN_VCF_IN_RESULTS            } from '../../../modules/local/replace_vcf_in_results/main'
include { REPLACE_VCF_IN_RESULTS as REPLACE_MSI_TN_VCF_IN_RESULTS            } from '../../../modules/local/replace_vcf_in_results/main'

include { BCFTOOLS_CONCAT                                                    } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_RENAME_TUMOR_SAMPLE         } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_SVC_TN_TO_ORIGINAL           } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_CNC_TN_TO_ORIGINAL           } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_MSI_TN_TO_ORIGINAL           } from '../../../modules/nf-core/bcftools/reheader/main'
include { XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS                           } from '../../../modules/local/xoos_copy_number_caller/tumor_normal_wgs/main'
include { XOOS_MSI_CALLER_TUMOR_NORMAL                                       } from '../../../modules/local/xoos_msi_caller/tumor_normal/main'

def removeBamSampleNames(meta) {
    def clean_meta = meta.clone()
    clean_meta.remove('tumor_sample_name')
    clean_meta.remove('normal_sample_name')
    clean_meta
}

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
    ch_pon_vcf_with_tbi
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
    ch_msi_sites_thresholds_tsv
    enable_msi_caller

    main:

    ch_all_vcfs = channel.empty()
    ch_raw_vc_vcf = channel.empty()
    ch_vcfeval_report = channel.empty()
    ch_vcfeval_output = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_copy_number_caller_output = channel.empty()
    ch_msi_caller_output = channel.empty()
    ch_tumor_normal_filtered_vcf = channel.empty()
    ch_mutect_stats = channel.empty()
    ch_tumor_normal_bams_with_sample_names = channel.empty()
    ch_tumor_sample_names = channel.empty()
    ch_svc_metrics = channel.empty()

    if (enable_small_variant_caller || enable_copy_number_caller || enable_msi_caller) {
        RESOLVE_TUMOR_NORMAL_BAM_SAMPLE_NAMES(ch_tumor_normal_bams)
        ch_tumor_normal_bams_with_sample_names = RESOLVE_TUMOR_NORMAL_BAM_SAMPLE_NAMES.out.bams
        ch_tumor_sample_names = ch_tumor_normal_bams_with_sample_names
            .map { meta, _bams, _bais -> [meta.id, meta.tumor_sample_name, meta.tumor_id ?: meta.id, meta.normal_sample_name, meta.normal_id ?: meta.id] }
    }

    // Step 1: Call small variants (Mutect2 + pretrained filtering)
    // Note: Somatic TN copy number calling requires the small variant VCF as input
    if (enable_small_variant_caller || enable_copy_number_caller) {
        CALL_SOMATIC_SMALL_VARIANTS_TN(
            ch_tumor_normal_bams_with_sample_names,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_pop_af_vcf_with_tbi,
            ch_svc_model_tumor_normal,
            ch_pon_vcf_with_tbi,
            ch_intervals,
            use_gatk_mutect2,
            enable_mutect_stats,
            enable_partitioning,
            partition_count,
            aligner,
        )
        ch_raw_vc_vcf = CALL_SOMATIC_SMALL_VARIANTS_TN.out.raw_vcf
        ch_tumor_normal_vcf_with_sample_names = CALL_SOMATIC_SMALL_VARIANTS_TN.out.vcf
            .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
            .join(ch_tumor_sample_names)
            .map { _sample_key, meta, vcf, tbi, tumor_sample_name, tumor_id, normal_sample_name, normal_id ->
                [meta + [tumor_sample_name: tumor_sample_name, tumor_id: tumor_id, normal_sample_name: normal_sample_name, normal_id: normal_id], vcf, tbi]
            }
        ch_svc_sample_names = ch_tumor_normal_vcf_with_sample_names.collectFile { meta, _vcf, _tbi ->
            ["${meta.id}.txt", "${meta.tumor_sample_name} ${meta.tumor_id}\n${meta.normal_sample_name} ${meta.normal_id}\n"]
        }
        ch_svc_reheader_input = ch_tumor_normal_vcf_with_sample_names
            .combine(ch_svc_sample_names)
            .filter { meta, _vcf, _tbi, sample_file -> sample_file.baseName == meta.id }
            .map { meta, vcf, _tbi, sample_file -> [meta, vcf, [], sample_file] }
        BCFTOOLS_REHEADER_SVC_TN_TO_ORIGINAL(
            ch_svc_reheader_input,
            channel.value([[], []]),
        )
        ch_tumor_normal_filtered_vcf = BCFTOOLS_REHEADER_SVC_TN_TO_ORIGINAL.out.vcf
            .join(BCFTOOLS_REHEADER_SVC_TN_TO_ORIGINAL.out.index)
            .map { meta, vcf, tbi -> [removeBamSampleNames(meta), vcf, tbi] }
        ch_mutect_stats = CALL_SOMATIC_SMALL_VARIANTS_TN.out.stats
        ch_svc_metrics = CALL_SOMATIC_SMALL_VARIANTS_TN.out.metrics
            .map { meta, metrics -> [meta + [svc_source: 'small_variant_tumor_normal'], metrics] }

        // Add pretrained VCFs to merged output
        ch_all_vcfs = ch_all_vcfs.mix(ch_tumor_normal_filtered_vcf)

        // Step 2: Evaluate somatic small variant calls with RTG vcfeval, if enabled
        if (enable_vcfeval) {
            if (enable_small_variant_caller) {
                // Filter out normal sample from the VCF input to vcfeval
                // and rename the tumor sample to match the truth VCF sample name
                // This is because vcfeval expects a single sample column with the same name
                // in both the query and truth VCFs.
                ch_vcfeval_candidates = ch_tumor_normal_vcf_with_sample_names.map { meta, vcf, tbi ->
                    [meta + [truth_sample_name: resolveTruthSampleName(meta, meta.id)], vcf, tbi]
                }
                ch_sample_names = ch_vcfeval_candidates.collectFile { meta, _vcf, _tbi ->
                    ["${meta.id}.txt", "${meta.tumor_sample_name} ${meta.truth_sample_name}" + '\n']
                }
                ch_bcftools_reheader_input = ch_vcfeval_candidates
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
                ).report.map { meta, report -> [removeBamSampleNames(meta), report] }.concat(ch_vcfeval_report)
                ch_multiqc_files = EVALUATE_SOMATIC_SMALL_VARIANTS_TN.out.custom_multiqc_files
                    .map { meta, multiqc_dir -> [removeBamSampleNames(meta), multiqc_dir] }
                    .concat(ch_multiqc_files)
                ch_vcfeval_output = EVALUATE_SOMATIC_SMALL_VARIANTS_TN.out.vcfeval_output
                    .map { meta, vcfeval_dir -> [removeBamSampleNames(meta), vcfeval_dir] }
                    .concat(ch_vcfeval_output)
            }
        }

        // Step 3: Call copy number variants with XOOS Copy Number Caller, if enabled
        if (enable_copy_number_caller) {
            ch_cnc_input = ch_tumor_normal_bams_with_sample_names
                .map { meta, tumor_normal_bams, tumor_normal_bam_indices -> [meta, tumor_normal_bams[0], tumor_normal_bam_indices[0], tumor_normal_bams[1], tumor_normal_bam_indices[1]] }
                .join(ch_tumor_normal_vcf_with_sample_names)
            XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS(
                ch_cnc_input,
                ch_reference,
                ch_mappability_bw,
                ch_cnv_seed_file,
                ch_cnv_blocklist_bed,
            )
            ch_copy_number_vcf_with_sample_names = XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS.out.vcf
            ch_cnc_sample_names = ch_copy_number_vcf_with_sample_names.collectFile { meta, _vcf ->
                ["${meta.id}.txt", "${meta.tumor_sample_name} ${meta.tumor_id}\n${meta.normal_sample_name} ${meta.normal_id}\n"]
            }
            ch_cnc_reheader_input = ch_copy_number_vcf_with_sample_names
                .combine(ch_cnc_sample_names)
                .filter { meta, _vcf, sample_file -> sample_file.baseName == meta.id }
                .map { meta, vcf, sample_file -> [meta, vcf, [], sample_file] }
            BCFTOOLS_REHEADER_CNC_TN_TO_ORIGINAL(
                ch_cnc_reheader_input,
                channel.value([[], []]),
            )
            ch_copy_number_vcf_with_original_names = BCFTOOLS_REHEADER_CNC_TN_TO_ORIGINAL.out.vcf.join(BCFTOOLS_REHEADER_CNC_TN_TO_ORIGINAL.out.index)
            ch_copy_number_vcf = ch_copy_number_vcf_with_original_names
                .map { meta, vcf, tbi -> [removeBamSampleNames(meta), vcf, tbi] }
            ch_copy_number_caller_output_input = XOOS_COPY_NUMBER_CALLER_TUMOR_NORMAL_WGS.out.results_dir
                .map { meta, results_dir -> [meta.id, meta, results_dir] }
                .join(ch_copy_number_vcf_with_original_names.map { meta, vcf, tbi -> [meta.id, vcf, tbi] })
                .map { _sample_key, meta, results_dir, vcf, tbi -> [meta, results_dir, vcf, tbi, '*_somatic_copy_number_callset.vcf.gz', 'copy_number_caller'] }
            REPLACE_CNC_TN_VCF_IN_RESULTS(ch_copy_number_caller_output_input)
            ch_copy_number_caller_output = REPLACE_CNC_TN_VCF_IN_RESULTS.out.results_dir
                .map { meta, results_dir -> [removeBamSampleNames(meta), results_dir] }
            ch_all_vcfs = ch_copy_number_vcf.mix(ch_all_vcfs)
        }
    }

    // Step 4: Call microsatellite instability (MSI) status, if enabled (beta).
    // MSI runs directly on the tumor/normal BAM pair and does not depend on the
    // small variant caller, so it is gated independently.
    if (enable_msi_caller) {
        ch_msi_input = ch_tumor_normal_bams_with_sample_names
            .map { meta, tumor_normal_bams, tumor_normal_bam_indices -> [meta, tumor_normal_bams[0], tumor_normal_bam_indices[0], tumor_normal_bams[1], tumor_normal_bam_indices[1]] }
        XOOS_MSI_CALLER_TUMOR_NORMAL(
            ch_msi_input,
            ch_msi_sites_thresholds_tsv,
        )
        ch_msi_vcf_with_sample_names = XOOS_MSI_CALLER_TUMOR_NORMAL.out.vcf
        ch_msi_sample_names = ch_msi_vcf_with_sample_names.collectFile { meta, _vcf ->
            ["${meta.id}.txt", "${meta.tumor_sample_name} ${meta.tumor_id}\n${meta.normal_sample_name} ${meta.normal_id}\n"]
        }
        ch_msi_reheader_input = ch_msi_vcf_with_sample_names
            .combine(ch_msi_sample_names)
            .filter { meta, _vcf, sample_file -> sample_file.baseName == meta.id }
            .map { meta, vcf, sample_file -> [meta, vcf, [], sample_file] }
        BCFTOOLS_REHEADER_MSI_TN_TO_ORIGINAL(
            ch_msi_reheader_input,
            channel.value([[], []]),
        )
        ch_msi_vcf_with_original_names = BCFTOOLS_REHEADER_MSI_TN_TO_ORIGINAL.out.vcf.join(BCFTOOLS_REHEADER_MSI_TN_TO_ORIGINAL.out.index)
        ch_msi_caller_output_input = XOOS_MSI_CALLER_TUMOR_NORMAL.out.results_dir
            .map { meta, results_dir -> [meta.id, meta, results_dir] }
            .join(ch_msi_vcf_with_original_names.map { meta, vcf, tbi -> [meta.id, vcf, tbi] })
            .map { _sample_key, meta, results_dir, vcf, tbi -> [meta, results_dir, vcf, tbi, '*_msi_results.vcf.gz', 'msi_caller'] }
        REPLACE_MSI_TN_VCF_IN_RESULTS(ch_msi_caller_output_input)
        ch_msi_caller_output = REPLACE_MSI_TN_VCF_IN_RESULTS.out.results_dir
            .map { meta, results_dir -> [removeBamSampleNames(meta), results_dir] }
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
    msi_caller_output            = ch_msi_caller_output // channel: [ val(meta), path(msi_results_dir) ]
    custom_multiqc_files         = ch_multiqc_files
    mutect_stats                 = ch_mutect_stats // channel: [ val(meta), path(stats) ]
    svc_metrics                  = ch_svc_metrics // channel: [ val(meta), path(metrics_tsv) ]
}
