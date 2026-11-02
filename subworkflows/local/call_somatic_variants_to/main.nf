include { CALL_SOMATIC_SMALL_VARIANTS_TO                                     } from './call_somatic_small_variants_to'
include { EVALUATE_SMALL_VARIANT_CALLS as EVALUATE_SOMATIC_SMALL_VARIANTS_TO } from '../evaluate_small_variant_calls/main'
include { RESOLVE_TUMOR_ONLY_BAM_SAMPLE_NAME                                } from './resolve_tumor_only_bam_sample_name'
include { REPLACE_VCF_IN_RESULTS as REPLACE_CNC_VCF_IN_RESULTS              } from '../../../modules/local/replace_vcf_in_results/main'
include { REPLACE_VCF_IN_RESULTS as REPLACE_MSI_VCF_IN_RESULTS              } from '../../../modules/local/replace_vcf_in_results/main'

include { BCFTOOLS_CONCAT                          } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_SVC_TO_ORIGINAL } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_CNC_TO_ORIGINAL } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_MSI_TO_ORIGINAL } from '../../../modules/nf-core/bcftools/reheader/main'
include { XOOS_COPY_NUMBER_CALLER_TUMOR_ONLY_WGS   } from '../../../modules/local/xoos_copy_number_caller/tumor_only_wgs/main'
include { XOOS_MSI_CALLER_TUMOR_ONLY               } from '../../../modules/local/xoos_msi_caller/tumor_only/main'
include { XOOS_GERMLINE_FILTERING                  } from '../../../modules/local/xoos_germline_filtering/main'

def removeBamSampleNames(meta) {
    def clean_meta = meta.clone()
    clean_meta.remove('tumor_sample_name')
    clean_meta.remove('normal_sample_name')
    clean_meta
}

workflow CALL_SOMATIC_VARIANTS_TO {
    take:
    ch_tumor_bams
    ch_reference
    ch_reference_dict
    ch_splitintervals_reference
    ch_splitintervals_reference_dict
    ch_intervals
    ch_pop_af_vcf_with_tbi
    ch_svc_model_tumor_only
    ch_pon_vcf_with_tbi
    ch_truth_vcf_dir_tumor_only
    ch_eval_bed_dir_tumor_only
    ch_mappability_bw
    ch_cnv_seed_file
    ch_cnv_blocklist_bed
    ch_cnv_snp_blocklist_bed
    enable_svc
    use_gatk_mutect2
    enable_mutect_stats
    enable_vcfeval
    enable_cnc
    enable_germline_filtering
    enable_merge_vcf
    reference_sdf
    enable_partitioning
    partition_count
    ch_msi_background_threshold_json
    enable_msi_caller

    main:

    ch_all_vcfs = channel.empty()
    ch_raw_vc_vcf = channel.empty()
    ch_vcfeval_report = channel.empty()
    ch_vcfeval_output = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_copy_number_caller_output = channel.empty()
    ch_tumor_only_filtered_vcf = channel.empty()
    ch_mutect_stats = channel.empty()
    ch_svc_metrics = channel.empty()
    ch_germline_filtering_output = channel.empty()
    ch_germline_filtering_vcf = channel.empty()
    ch_msi_caller_output = channel.empty()
    ch_tumor_bams_with_sample_names = channel.empty()
    ch_tumor_sample_names = channel.empty()

    if (enable_svc || enable_cnc || enable_msi_caller) {
        RESOLVE_TUMOR_ONLY_BAM_SAMPLE_NAME(ch_tumor_bams)
        ch_tumor_bams_with_sample_names = RESOLVE_TUMOR_ONLY_BAM_SAMPLE_NAME.out.bams
        ch_tumor_sample_names = ch_tumor_bams_with_sample_names
            .map { meta, _bams, _bais -> [meta.id, meta.tumor_sample_name, meta.tumor_id ?: meta.id] }
    }

    // Step 1: Call small variants (Mutect2 tumor-only + pretrained filtering)
    // Note: Somatic TO copy number calling requires the small variant VCF as input
    if (enable_svc || enable_cnc) {
        CALL_SOMATIC_SMALL_VARIANTS_TO(
            ch_tumor_bams_with_sample_names,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_pop_af_vcf_with_tbi,
            ch_svc_model_tumor_only,
            ch_pon_vcf_with_tbi,
            ch_intervals,
            use_gatk_mutect2,
            enable_mutect_stats,
            enable_partitioning,
            partition_count,
        )
        ch_raw_vc_vcf = CALL_SOMATIC_SMALL_VARIANTS_TO.out.raw_vcf
        ch_mutect_stats = CALL_SOMATIC_SMALL_VARIANTS_TO.out.stats
        ch_svc_metrics = CALL_SOMATIC_SMALL_VARIANTS_TO.out.metrics
            .map { meta, metrics -> [meta + [svc_source: 'small_variant_tumor_only'], metrics] }

        // Reattach the resolved SM metadata to the public SVC output from the
        // nested workflow before using it in downstream commands.
        ch_svc_vcf_with_sample_names = CALL_SOMATIC_SMALL_VARIANTS_TO.out.vcf
            .map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] }
            .join(ch_tumor_sample_names)
            .map { _sample_key, meta, vcf, tbi, tumor_sample_name, tumor_id ->
                [meta + [tumor_sample_name: tumor_sample_name, tumor_id: tumor_id], vcf, tbi]
            }

        // Reheader the final SVC VCF back to the sample-sheet tumor name. The
        // SM value remains in the internal channel used by RTG and CNC.
        ch_svc_sample_names = ch_svc_vcf_with_sample_names.collectFile { meta, _vcf, _tbi ->
            ["${meta.id}.txt", "${meta.tumor_sample_name} ${meta.tumor_id}" + '\n']
        }
        ch_reheader_input = ch_svc_vcf_with_sample_names
            .combine(ch_svc_sample_names)
            .filter { meta, _vcf, _tbi, sample_file -> sample_file.baseName == meta.id }
            .map { meta, vcf, _tbi, sample_file -> [meta, vcf, [], sample_file] }
        BCFTOOLS_REHEADER_SVC_TO_ORIGINAL(
            ch_reheader_input,
            channel.value([[], []]),
        )
        ch_tumor_only_filtered_vcf_with_original_names = BCFTOOLS_REHEADER_SVC_TO_ORIGINAL.out.vcf.join(BCFTOOLS_REHEADER_SVC_TO_ORIGINAL.out.index)
        ch_tumor_only_filtered_vcf_with_sample_names = ch_svc_vcf_with_sample_names
        ch_tumor_only_filtered_vcf = ch_tumor_only_filtered_vcf_with_original_names
            .map { meta, vcf, tbi -> [removeBamSampleNames(meta), vcf, tbi] }

        // Step 2: Call copy number variants with XOOS Copy Number Caller, if enabled
        if (enable_cnc) {
            ch_cnc_input = ch_tumor_bams_with_sample_names
                // The resolver adds the BAM SM tag as `tumor_sample_name`, while
                // the reattached SVC metadata also adds the fallback `tumor_id`.
                // Join by the stable sample id rather than the complete metadata
                // map, otherwise samples whose input metadata lacks tumor_id are
                // silently dropped before CNC is scheduled.
                .map { meta, tumor_bams, tumor_bam_indices -> [meta.id, tumor_bams[0], tumor_bam_indices[0]] }
                .join(ch_tumor_only_filtered_vcf_with_sample_names.map { meta, vcf, tbi -> [meta.id, meta, vcf, tbi] })
                .map { _sample_key, tumor_bam, tumor_bam_index, meta, vcf, tbi -> [meta, tumor_bam, tumor_bam_index, vcf, tbi] }
            XOOS_COPY_NUMBER_CALLER_TUMOR_ONLY_WGS(
                ch_cnc_input,
                ch_reference,
                ch_mappability_bw,
                ch_cnv_seed_file,
                ch_cnv_blocklist_bed,
                ch_pop_af_vcf_with_tbi,
                ch_cnv_snp_blocklist_bed,
            )
            ch_copy_number_vcf_with_sample_names = XOOS_COPY_NUMBER_CALLER_TUMOR_ONLY_WGS.out.vcf
            ch_cnc_sample_names = ch_copy_number_vcf_with_sample_names.collectFile { meta, _vcf ->
                ["${meta.id}.txt", "${meta.tumor_sample_name} ${meta.tumor_id}" + '\n']
            }
            ch_cnc_reheader_input = ch_copy_number_vcf_with_sample_names
                .combine(ch_cnc_sample_names)
                .filter { meta, _vcf, sample_file -> sample_file.baseName == meta.id }
                .map { meta, vcf, sample_file -> [meta, vcf, [], sample_file] }
            BCFTOOLS_REHEADER_CNC_TO_ORIGINAL(
                ch_cnc_reheader_input,
                channel.value([[], []]),
            )
            ch_copy_number_vcf_with_original_names = BCFTOOLS_REHEADER_CNC_TO_ORIGINAL.out.vcf.join(BCFTOOLS_REHEADER_CNC_TO_ORIGINAL.out.index)
            ch_copy_number_vcf = ch_copy_number_vcf_with_original_names
                .map { meta, vcf, tbi -> [removeBamSampleNames(meta), vcf, tbi] }
            ch_copy_number_caller_output_input = XOOS_COPY_NUMBER_CALLER_TUMOR_ONLY_WGS.out.results_dir
                .map { meta, results_dir -> [meta.id, meta, results_dir] }
                .join(ch_copy_number_vcf_with_original_names.map { meta, vcf, tbi -> [meta.id, vcf, tbi] })
                .map { _sample_key, meta, results_dir, vcf, tbi -> [meta, results_dir, vcf, tbi, '*_somatic_copy_number_callset.vcf.gz', 'copy_number_caller'] }
            REPLACE_CNC_VCF_IN_RESULTS(ch_copy_number_caller_output_input)
            ch_copy_number_caller_output = REPLACE_CNC_VCF_IN_RESULTS.out.results_dir
                .map { meta, results_dir -> [removeBamSampleNames(meta), results_dir] }
            ch_copy_number_seg = XOOS_COPY_NUMBER_CALLER_TUMOR_ONLY_WGS.out.seg
                .map { meta, seg -> [removeBamSampleNames(meta), seg] }
            ch_copy_number_vcf_with_index = ch_copy_number_vcf

            // Step 4: Germline filtering, consuming CNC purity/ploidy (SEG).
            // SNV/INDEL models default to the ones shipped in the container.
            if (enable_germline_filtering) {
                ch_germline_filtering_input = ch_tumor_only_filtered_vcf
                    .map { meta, vcf, _tbi -> [meta, vcf] }
                    .join(ch_copy_number_seg)
                XOOS_GERMLINE_FILTERING(
                    ch_germline_filtering_input,
                )
                ch_germline_filtering_output = XOOS_GERMLINE_FILTERING.out.prediction
                ch_germline_filtering_vcf = XOOS_GERMLINE_FILTERING.out.vcf
            }
        }

        // Step 3: Evaluate somatic small variant calls with RTG vcfeval, if enabled.
        // When germline filtering is enabled, evaluate its tagged VCF because it
        // is the small-variant output used for the tumor-only workflow. Otherwise,
        // evaluate the SVC filter-variants VCF directly.
        if (enable_vcfeval && enable_svc) {
            ch_vcfeval_input = enable_germline_filtering
                ? ch_germline_filtering_vcf
                : ch_tumor_only_filtered_vcf
            ch_vcfeval_report = EVALUATE_SOMATIC_SMALL_VARIANTS_TO(
                ch_vcfeval_input.map { meta, vcf, tbi -> [meta + [svc_model_type: 'small_variant_pretrained'], vcf, tbi] },
                ch_reference,
                ch_truth_vcf_dir_tumor_only,
                ch_eval_bed_dir_tumor_only,
                reference_sdf,
            ).report.map { meta, report -> [removeBamSampleNames(meta), report] }.concat(ch_vcfeval_report)
            ch_multiqc_files = EVALUATE_SOMATIC_SMALL_VARIANTS_TO.out.custom_multiqc_files
                .map { meta, multiqc_dir -> [removeBamSampleNames(meta), multiqc_dir] }
                .concat(ch_multiqc_files)
            ch_vcfeval_output = EVALUATE_SOMATIC_SMALL_VARIANTS_TO.out.vcfeval_output
                .map { meta, vcfeval_dir -> [removeBamSampleNames(meta), vcfeval_dir] }
                .concat(ch_vcfeval_output)
        }

        // Assemble the small-variant contribution to the merge. When germline
        // filtering runs, its VCF supersedes the raw small-variant VCF;
        // otherwise the small-variant VCF is used directly.
        ch_small_variant_for_merge = enable_germline_filtering
            ? ch_germline_filtering_vcf
            : ch_tumor_only_filtered_vcf
        ch_all_vcfs = ch_all_vcfs.mix(ch_small_variant_for_merge)

        if (enable_cnc) {
            ch_all_vcfs = ch_all_vcfs.mix(ch_copy_number_vcf_with_index)
        }
    }

    if (enable_msi_caller) {
        ch_msi_input = ch_tumor_bams_with_sample_names
            .map { meta, tumor_bams, tumor_bam_indices -> [meta, tumor_bams[0], tumor_bam_indices[0]] }
        XOOS_MSI_CALLER_TUMOR_ONLY(ch_msi_input, ch_msi_background_threshold_json)
        ch_msi_vcf_with_sample_names = XOOS_MSI_CALLER_TUMOR_ONLY.out.vcf
            // Keep the MSI metadata aligned with the SVC/CNC contribution when
            // the input samplesheet does not provide tumor_id. The SVC path
            // uses the sample id as the fallback, so the merge group must use
            // the same canonical metadata for all caller VCFs.
            .map { meta, vcf -> [meta + [tumor_id: meta.tumor_id ?: meta.id], vcf] }
        ch_msi_sample_names = ch_msi_vcf_with_sample_names.collectFile { meta, _vcf ->
            ["${meta.id}.txt", "${meta.tumor_sample_name} ${meta.tumor_id}" + '\n']
        }
        ch_msi_reheader_input = ch_msi_vcf_with_sample_names
            .combine(ch_msi_sample_names)
            .filter { meta, _vcf, sample_file -> sample_file.baseName == meta.id }
            .map { meta, vcf, sample_file -> [meta, vcf, [], sample_file] }
        BCFTOOLS_REHEADER_MSI_TO_ORIGINAL(
            ch_msi_reheader_input,
            channel.value([[], []]),
        )
        ch_msi_vcf_with_original_names = BCFTOOLS_REHEADER_MSI_TO_ORIGINAL.out.vcf.join(BCFTOOLS_REHEADER_MSI_TO_ORIGINAL.out.index)
        ch_msi_caller_output_input = XOOS_MSI_CALLER_TUMOR_ONLY.out.results_dir
            .map { meta, results_dir -> [meta.id, meta, results_dir] }
            .join(ch_msi_vcf_with_original_names.map { meta, vcf, tbi -> [meta.id, vcf, tbi] })
            .map { _sample_key, meta, results_dir, vcf, tbi -> [meta, results_dir, vcf, tbi, '*_msi_results.vcf.gz', 'msi_caller'] }
        REPLACE_MSI_VCF_IN_RESULTS(ch_msi_caller_output_input)
        ch_msi_caller_output = REPLACE_MSI_VCF_IN_RESULTS.out.results_dir
            .map { meta, results_dir -> [removeBamSampleNames(meta), results_dir] }

        // Include the reheadered MSI calls in the merged tumor-only VCF. The
        // MSI results directory is published separately above, but omitting
        // this VCF from ch_all_vcfs makes enable_merge_vcf silently drop MSI
        // records from the sample-level output.
        ch_all_vcfs = ch_all_vcfs.mix(
            ch_msi_vcf_with_original_names.map { meta, vcf, tbi ->
                [removeBamSampleNames(meta), vcf, tbi]
            }
        )
    }

    ch_merged_vcf = channel.empty()
    if (enable_merge_vcf) {
        // Count enabled callers so groupTuple knows the expected group size per sample.
        // SVC is run whenever CNC is enabled, while MSI is independently
        // optional. Count only the VCF-producing callers that feed the merge.
        def num_callers = (enable_svc || enable_cnc ? 1 : 0) +
            (enable_cnc ? 1 : 0) +
            (enable_msi_caller ? 1 : 0)
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
    svc_vcf                      = ch_tumor_only_filtered_vcf // [ meta, path(vcf), path(tbi) ]
    // The germline-filtered VCF (non-empty only when enable_germline_filtering is true).
    // Exposed separately so callers can publish it in the no-merge path.
    germline_filtered_vcf        = ch_germline_filtering_vcf // [ meta, path(vcf), path(tbi) ]
    small_variant_vcfeval_report = ch_vcfeval_report
    small_variant_vcfeval_output = ch_vcfeval_output // channel: [ val(meta), path(vcfeval_dir) ]
    copy_number_caller_output    = ch_copy_number_caller_output
    germline_filtering_output    = ch_germline_filtering_output // channel: [ val(meta), path(prediction) ]
    msi_caller_output            = ch_msi_caller_output // channel: [ val(meta), path(msi_results_dir) ]
    custom_multiqc_files         = ch_multiqc_files
    mutect_stats                 = ch_mutect_stats // channel: [ val(meta), path(stats) ]
    svc_metrics                  = ch_svc_metrics // channel: [ val(meta), path(metrics_tsv) ]
}
