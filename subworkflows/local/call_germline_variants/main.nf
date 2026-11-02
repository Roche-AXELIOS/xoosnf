include { CALL_GERMLINE_SMALL_VARIANTS                                                } from './call_germline_small_variants'
include { EVALUATE_SMALL_VARIANT_CALLS as EVALUATE_GERMLINE_SMALL_VARIANTS_PRETRAINED ; EVALUATE_SMALL_VARIANT_CALLS as EVALUATE_GERMLINE_SMALL_VARIANTS_RETRAINED } from '../evaluate_small_variant_calls/main'

include { XOOS_COPY_NUMBER_CALLER_GERMLINE_WGS                                        } from '../../../modules/local/xoos_copy_number_caller/germline_wgs/main'
include { RESOLVE_GERMLINE_BAM_SAMPLE_NAME                                       } from './resolve_germline_bam_sample_name'
include { REPLACE_VCF_IN_RESULTS as REPLACE_CNC_GERMLINE_VCF_IN_RESULTS          } from '../../../modules/local/replace_vcf_in_results/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_GERMLINE_SAMPLE ; BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_GERMLINE_RETRAINED_SAMPLE } from '../../../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_CNC_GERMLINE                    } from '../../../modules/nf-core/bcftools/reheader/main'
include { XOOS_STR_CALLER                                                             } from '../../../modules/local/xoos_str_caller/main'
include { MITOCHONDRIAL_VARIANT_CALLING                                                } from '../mitochondrial_variant_calling/main'
include { BCFTOOLS_CONCAT                                                             } from '../../../modules/nf-core/bcftools/concat'

def removeBamSampleName(meta) {
    def clean_meta = meta.clone()
    clean_meta.remove('bam_sample_name')
    clean_meta.remove('germline_sample_name')
    clean_meta
}

workflow CALL_GERMLINE_VARIANTS {
    take:
    ch_bam
    ch_reference
    ch_reference_dict
    ch_splitintervals_reference
    ch_splitintervals_reference_dict
    ch_intervals
    ch_par_beds
    ch_pop_af_vcf_with_tbi
    ch_svc_models
    ch_truth_vcf_dir
    ch_train_bed_dir
    ch_test_bed_dir
    ch_eval_bed_dir
    ch_mappability_bw
    ch_cnv_seed_file
    ch_cnv_blocklist_bed
    ch_repeat_catalog_json
    enable_small_variant_caller
    use_gatk_haplotype_caller
    aligner
    enable_small_variant_caller_retrain
    enable_vcfeval
    enable_copy_number_caller
    enable_copy_number_caller_baf
    enable_str_caller
    enable_mitochondrial_variant_caller
    use_gatk_mutect2_mitochondrial
    enable_merge_vcf
    reference_sdf
    enable_partitioning
    partition_count
    ch_giraffe_graphs
    ch_giraffe_chrm_paths

    main:

    ch_all_vcfs = channel.empty()
    ch_raw_vc_vcf = channel.empty()
    ch_small_variant_vcf = channel.empty()
    ch_small_variant_retrained_vcf = channel.empty()
    ch_copy_number_caller_output = channel.empty()
    ch_str_caller_output = channel.empty()
    ch_mitochondrial_vcf = channel.empty()
    ch_vcfeval_report = channel.empty()
    ch_vcfeval_output = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_svc_metrics = channel.empty()
    ch_bams_with_sample_names = channel.empty()

    // The BAM SM tag is required only when CNC BAF is enabled, because that
    // mode passes the HaplotypeCaller VCF into CNC and therefore needs CNC's
    // internal sample ID to match the BAM-derived VCF sample column. For
    // ordinary CNC calls, the samplesheet ID remains the sample name and no
    // BAM header inspection is needed.
    if (enable_copy_number_caller && enable_copy_number_caller_baf) {
        RESOLVE_GERMLINE_BAM_SAMPLE_NAME(ch_bam)
        ch_bams_with_sample_names = RESOLVE_GERMLINE_BAM_SAMPLE_NAME.out.bams
    } else {
        ch_bams_with_sample_names = ch_bam
    }

    // Step 1: Call small variants (HaplotypeCaller + pretrained/retrained filtering).
    // Also runs when CNC BAF is requested so HaplotypeCaller VCF is available for BAF.
    if (enable_small_variant_caller || enable_small_variant_caller_retrain || (enable_copy_number_caller && enable_copy_number_caller_baf)) {
        CALL_GERMLINE_SMALL_VARIANTS(
            ch_bam,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_par_beds,
            ch_pop_af_vcf_with_tbi,
            ch_svc_models,
            ch_truth_vcf_dir,
            ch_train_bed_dir,
            ch_test_bed_dir,
            ch_intervals,
            use_gatk_haplotype_caller,
            aligner,
            enable_small_variant_caller,
            enable_small_variant_caller_retrain,
            enable_partitioning,
            partition_count,
        )

        ch_raw_vc_vcf = CALL_GERMLINE_SMALL_VARIANTS.out.raw_vcf
        ch_svc_metrics = CALL_GERMLINE_SMALL_VARIANTS.out.metrics_pretrained_model
            .map { meta, metrics -> [meta + [svc_source: 'small_variant_pretrained'], metrics] }
            .mix(CALL_GERMLINE_SMALL_VARIANTS.out.metrics_retrained_model
                .map { meta, metrics -> [meta + [svc_source: 'small_variant_retrained'], metrics] })

        // Add pretrained VCFs to merged output
        ch_pretrained_filtered_vcfs = channel.empty()
        if (enable_small_variant_caller) {
            // Reheader the SVC VCF to set the sample name to meta.id.
            // BAMs often carry an incorrect SM read-group tag, which the small
            // variant caller propagates into the VCF sample column. Reheadering
            // here ensures the output VCF always uses the canonical sample name.
            ch_svc_sample_names = CALL_GERMLINE_SMALL_VARIANTS.out.vcf_pretrained_model.collectFile { meta, _vcf, _tbi ->
                ["${meta.id}.txt", "${meta.id}" + '\n']
            }
            ch_bcftools_reheader_input = CALL_GERMLINE_SMALL_VARIANTS.out.vcf_pretrained_model
                .combine(ch_svc_sample_names)
                .filter { meta, _vcf, _tbi, sample_file -> sample_file.baseName == meta.id }
                .map { meta, vcf, _tbi, sample_file -> [meta, vcf, [], sample_file] }
            BCFTOOLS_REHEADER_GERMLINE_SAMPLE(
                ch_bcftools_reheader_input,
                channel.value([[], []]),
            )
            ch_small_variant_vcf = BCFTOOLS_REHEADER_GERMLINE_SAMPLE.out.vcf.join(BCFTOOLS_REHEADER_GERMLINE_SAMPLE.out.index)

            ch_pretrained_filtered_vcfs = ch_small_variant_vcf
                .map { meta, vcf, tbi -> [meta, 'small_variant_pretrained', vcf, tbi] }
                .concat(ch_pretrained_filtered_vcfs)
            ch_all_vcfs = ch_pretrained_filtered_vcfs.concat(ch_all_vcfs)
        }

        // Add retrained VCFs to merged output
        ch_retrained_vcfs = channel.empty()
        if (enable_small_variant_caller_retrain) {
            // Reheader retrained VCF for the same reason as the pretrained VCF above.
            ch_retrained_sample_names = CALL_GERMLINE_SMALL_VARIANTS.out.vcf_retrained_model.collectFile { meta, _vcf, _tbi ->
                ["${meta.id}.txt", "${meta.id}" + '\n']
            }
            ch_retrained_reheader_input = CALL_GERMLINE_SMALL_VARIANTS.out.vcf_retrained_model
                .combine(ch_retrained_sample_names)
                .filter { meta, _vcf, _tbi, sample_file -> sample_file.baseName == meta.id }
                .map { meta, vcf, _tbi, sample_file -> [meta, vcf, [], sample_file] }
            BCFTOOLS_REHEADER_GERMLINE_RETRAINED_SAMPLE(
                ch_retrained_reheader_input,
                channel.value([[], []]),
            )
            ch_retrained_reheadered = BCFTOOLS_REHEADER_GERMLINE_RETRAINED_SAMPLE.out.vcf.join(BCFTOOLS_REHEADER_GERMLINE_RETRAINED_SAMPLE.out.index)
            ch_small_variant_retrained_vcf = ch_retrained_reheadered

            ch_retrained_vcfs = ch_retrained_reheadered
                .map { meta, vcf, tbi -> [meta, 'small_variant_retrained', vcf, tbi] }
                .concat(ch_retrained_vcfs)
        }

        // Evaluate germline small variant calls with RTG vcfeval, if enabled
        if (enable_vcfeval) {
            if (enable_small_variant_caller) {
                ch_vcfeval_report = EVALUATE_GERMLINE_SMALL_VARIANTS_PRETRAINED(
                    ch_pretrained_filtered_vcfs.map { meta, type, vcf, tbi -> [meta + [svc_model_type: type], vcf, tbi] },
                    ch_reference,
                    ch_truth_vcf_dir,
                    ch_eval_bed_dir,
                    reference_sdf,
                ).report.concat(ch_vcfeval_report)
                ch_multiqc_files = EVALUATE_GERMLINE_SMALL_VARIANTS_PRETRAINED.out.custom_multiqc_files.concat(ch_multiqc_files)
                ch_vcfeval_output = EVALUATE_GERMLINE_SMALL_VARIANTS_PRETRAINED.out.vcfeval_output.concat(ch_vcfeval_output)
            }
            if (enable_small_variant_caller_retrain) {
                // For retrained model, use the test BEDs as evaluation regions
                // since they are trained on the train BED regions and we want an unbiased evaluation
                ch_vcfeval_report = EVALUATE_GERMLINE_SMALL_VARIANTS_RETRAINED(
                    ch_retrained_vcfs.map { meta, type, vcf, tbi -> [meta + [svc_model_type: type], vcf, tbi] },
                    ch_reference,
                    ch_truth_vcf_dir,
                    ch_test_bed_dir,
                    reference_sdf,
                ).report.concat(ch_vcfeval_report)
                ch_vcfeval_output = EVALUATE_GERMLINE_SMALL_VARIANTS_RETRAINED.out.vcfeval_output.concat(ch_vcfeval_output)
            }
        }
    }

    if (enable_copy_number_caller) {
        // CNC needs the BAM's SM tag to match the sample column in the raw
        // HaplotypeCaller VCF when BAF is enabled. The resolver retains the
        // samplesheet ID as the stable channel key for joins.
        ch_cnc_bams = ch_bams_with_sample_names
            .map { meta, bam, bai -> [meta.id, meta, bam, bai] }

        // Use the HaplotypeCaller raw VCF for BAF — not the SVC-filtered output —
        // so BAF data is available even when enable_small_variant_caller is false.
        // ch_raw_vc_vcf is guaranteed non-empty here because the
        // CALL_GERMLINE_SMALL_VARIANTS guard above runs whenever enable_copy_number_caller_baf is true.
        ch_cnc_input = enable_copy_number_caller_baf
            ? ch_cnc_bams
                // The resolver enriches BAM metadata with bam_sample_name, so
                // join on the stable samplesheet ID instead of the full meta map.
                .join(ch_raw_vc_vcf.map { meta, vcf, tbi -> [meta.id, vcf, tbi] }, by: 0)
                .map { _sample_key, meta, bam, bai, vcf, tbi -> [meta, bam, bai, vcf, tbi] }
            : ch_bams_with_sample_names.map { meta, bam, bai -> [meta, bam, bai, [], []] }
        XOOS_COPY_NUMBER_CALLER_GERMLINE_WGS(
            ch_cnc_input,
            ch_reference,
            ch_mappability_bw,
            ch_cnv_seed_file,
            ch_cnv_blocklist_bed,
        )
        ch_cnv_vcfs_with_bam_sample_names = XOOS_COPY_NUMBER_CALLER_GERMLINE_WGS.out.vcf
        ch_cnc_sample_names = ch_cnv_vcfs_with_bam_sample_names.collectFile { meta, _vcf ->
            ["${meta.id}.txt", "${meta.bam_sample_name ?: meta.germline_sample_name ?: meta.sample_name ?: meta.id} ${meta.id}" + '\n']
        }
        ch_cnc_reheader_input = ch_cnv_vcfs_with_bam_sample_names
            .combine(ch_cnc_sample_names)
            .filter { meta, _vcf, sample_file -> sample_file.baseName == meta.id }
            .map { meta, vcf, sample_file -> [meta, vcf, [], sample_file] }
        BCFTOOLS_REHEADER_CNC_GERMLINE(
            ch_cnc_reheader_input,
            channel.value([[], []]),
        )
        ch_cnv_vcfs_reheadered = BCFTOOLS_REHEADER_CNC_GERMLINE.out.vcf
            .join(BCFTOOLS_REHEADER_CNC_GERMLINE.out.index)

        ch_cnc_output_input = XOOS_COPY_NUMBER_CALLER_GERMLINE_WGS.out.results_dir
            .map { meta, results_dir -> [meta.id, meta, results_dir] }
            // Reheader output retains the enriched meta map; use the same
            // stable ID key when pairing it with the caller results directory.
            .join(ch_cnv_vcfs_reheadered.map { meta, vcf, tbi -> [meta.id, vcf, tbi] }, by: 0)
            .map { _sample_key, meta, results_dir, vcf, tbi ->
                [meta, results_dir, vcf, tbi, '*_germline_copy_number_callset.vcf.gz', 'copy_number_caller']
            }
        REPLACE_CNC_GERMLINE_VCF_IN_RESULTS(ch_cnc_output_input)
        ch_copy_number_caller_output = REPLACE_CNC_GERMLINE_VCF_IN_RESULTS.out.results_dir
            .map { meta, results_dir -> [removeBamSampleName(meta), results_dir] }
        ch_cnv_vcfs = ch_cnv_vcfs_reheadered
            .map { meta, vcf, tbi -> [removeBamSampleName(meta), vcf, tbi] }
        // Add the reheadered CNV VCF to merged output.
        ch_all_vcfs = ch_cnv_vcfs
            .map { meta, vcf, tbi -> [meta, 'copy_number', vcf, tbi] }
            .concat(ch_all_vcfs)
    }

    if (enable_str_caller) {
        XOOS_STR_CALLER(
            ch_bam,
            ch_reference,
            ch_repeat_catalog_json,
        )
        ch_str_vcfs = XOOS_STR_CALLER.out.vcf
        ch_str_caller_output = XOOS_STR_CALLER.out.results_dir
        // Add STR VCFs to merged output
        ch_all_vcfs = ch_str_vcfs
            .map { meta, vcf, tbi -> [meta, 'str', vcf, tbi] }
            .concat(ch_all_vcfs)
    }

    // Call mitochondrial variants (Mutect2 mito-mode + SVC germline-mitochondria)
    if (enable_mitochondrial_variant_caller) {
        MITOCHONDRIAL_VARIANT_CALLING(
            ch_bam,
            ch_reference,
            ch_reference_dict,
            use_gatk_mutect2_mitochondrial,
            ch_giraffe_graphs,
            ch_giraffe_chrm_paths,
        )
        ch_mitochondrial_vcf = MITOCHONDRIAL_VARIANT_CALLING.out.vcf
        ch_svc_metrics = ch_svc_metrics.mix(
            MITOCHONDRIAL_VARIANT_CALLING.out.metrics
                .map { meta, metrics -> [meta + [svc_source: 'small_variant_mitochondria'], metrics] }
        )
        // Add mitochondrial VCFs to merged output
        ch_all_vcfs = ch_mitochondrial_vcf
            .map { meta, vcf, tbi -> [meta, 'mitochondria', vcf, tbi] }
            .concat(ch_all_vcfs)
    }

    ch_merged_vcf = channel.empty()
    if (enable_merge_vcf) {
        // Count enabled callers so groupTuple knows the expected group size per sample.
        // Retrained SVC VCFs are intentionally excluded from the merge — they are
        // published separately for benchmarking and should not appear in the merged output.
        def num_callers = (enable_small_variant_caller ? 1 : 0) +
            (enable_copy_number_caller ? 1 : 0) +
            (enable_str_caller ? 1 : 0) +
            (enable_mitochondrial_variant_caller ? 1 : 0)
        ch_merged_vcf = BCFTOOLS_CONCAT(
            ch_all_vcfs
                .map { meta, type, vcf, tbi -> [groupKey(meta, num_callers), type, vcf, tbi] }
                .groupTuple()
                .map { meta, _type, vcf, tbi -> [meta, vcf, tbi] }
        ).vcf
        ch_merged_vcf = ch_merged_vcf.join(BCFTOOLS_CONCAT.out.index)
    } else {
        // Emit individual VCFs without merging
        ch_merged_vcf = ch_all_vcfs.map { meta, _type, vcf, tbi -> [meta, vcf, tbi] }
    }

    emit:
    // When enable_merge_vcf is true: one [meta, vcf, tbi] tuple per sample (merged).
    // When enable_merge_vcf is false: one [meta, vcf, tbi] tuple per caller per sample (unmerged).
    merged_vcf                   = ch_merged_vcf
    raw_vc_vcf                   = ch_raw_vc_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    svc_pretrained_vcf           = ch_small_variant_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    svc_retrained_vcf            = ch_small_variant_retrained_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    copy_number_caller_output    = ch_copy_number_caller_output
    str_caller_output            = ch_str_caller_output
    mitochondrial_vcf            = ch_mitochondrial_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    small_variant_vcfeval_report = ch_vcfeval_report
    small_variant_vcfeval_output = ch_vcfeval_output // channel: [ val(meta), path(vcfeval_dir) ]
    svc_metrics               = ch_svc_metrics // channel: [ val(meta), path(metrics_tsv) ]
    custom_multiqc_files         = ch_multiqc_files
}
