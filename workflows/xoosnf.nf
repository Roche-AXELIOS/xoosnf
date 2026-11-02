/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DEMULTIPLEX                    } from '../subworkflows/local/demultiplex/main'
include { ALIGN                          } from '../subworkflows/local/align/main'
include { PARTITION_ALIGN_MERGE          } from '../subworkflows/local/align/partition_align_merge.nf'
include { DEDUP                          } from '../subworkflows/local/dedup'
include { SUBSAMPLE_BAM                  } from '../subworkflows/local/subsample_bam'
include { ROUTE_INPUTS                   } from '../subworkflows/local/route_inputs/main'
include { METRICS as PRE_DEDUP_METRICS ; METRICS as DEDUP_METRICS ; METRICS as SUBSAMPLE_METRICS } from '../subworkflows/local/metrics/main'
include { CONTAMINATION_DETECTION        } from '../subworkflows/local/metrics/contamination_detection.nf'
include { TUMOR_FRACTION_ESTIMATION      } from '../subworkflows/local/tumor_fraction_estimation/main'
include { CALL_GERMLINE_VARIANTS         } from '../subworkflows/local/call_germline_variants/main'
include { CALL_SOMATIC_VARIANTS_TN       } from '../subworkflows/local/call_somatic_variants_tn/main'
include { SOMATIC_TN_PAIR_VERIFICATION   } from '../subworkflows/local/somatic_tn_pair_verification/main'
include { CALL_SOMATIC_VARIANTS_TO       } from '../subworkflows/local/call_somatic_variants_to/main'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { XOOS_GENERATE_REPORT           } from '../modules/local/xoos_generate_metrics_report/generate_report/main'
include { XOOS_GENERATE_REPORT_SUMMARY   } from '../modules/local/xoos_generate_metrics_report/generate_report_summary/main'
include { XOOS_PARSE_CONTAMINATION_METRICS                                   } from '../modules/local/xoos_generate_metrics_report/parse_contamination_metrics/main'
include { XOOS_PARSE_TUMOR_FRACTION_METRICS                                   } from '../modules/local/xoos_generate_metrics_report/parse_tumor_fraction_metrics/main'
include { XOOS_PARSE_STR_METRICS                                               } from '../modules/local/xoos_generate_metrics_report/parse_str_metrics/main'
include { XOOS_PARSE_SVC_METRICS                                               } from '../modules/local/xoos_generate_metrics_report/parse_svc_metrics/main'
include { XOOS_PARSE_CNC_METRICS as XOOS_PARSE_CNC_METRICS_SOMATIC            } from '../modules/local/xoos_generate_metrics_report/parse_cnc_metrics/main'
include { XOOS_PARSE_CNC_METRICS as XOOS_PARSE_CNC_METRICS_SOMATIC_TO         } from '../modules/local/xoos_generate_metrics_report/parse_cnc_metrics/main'
include { XOOS_PARSE_CNC_METRICS as XOOS_PARSE_CNC_METRICS_GERMLINE           } from '../modules/local/xoos_generate_metrics_report/parse_cnc_metrics/main'
include { XOOS_COMBINE_MULTIQC_REPORTS   } from '../modules/local/xoos_metrics_to_multiqc/combine_reports/main'
include { XOOS_STR_METRICS_TO_MULTIQC                                        } from '../modules/local/xoos_metrics_to_multiqc/str_metrics/main'
include { XOOS_CNC_METRICS_TO_MULTIQC as XOOS_CNC_METRICS_TO_MULTIQC_SOMATIC     } from '../modules/local/xoos_metrics_to_multiqc/cnc_metrics/main'
include { XOOS_CNC_METRICS_TO_MULTIQC as XOOS_CNC_METRICS_TO_MULTIQC_SOMATIC_TO } from '../modules/local/xoos_metrics_to_multiqc/cnc_metrics/main'
include { XOOS_CNC_METRICS_TO_MULTIQC as XOOS_CNC_METRICS_TO_MULTIQC_GERMLINE   } from '../modules/local/xoos_metrics_to_multiqc/cnc_metrics/main'
include { XOOS_READ_COLLAPSER_METRICS_TO_MULTIQC  } from '../modules/local/xoos_metrics_to_multiqc/read_collapser_metrics/main'
include { ADD_RUN_INFO_TO_MULTIQC_HEADER } from '../modules/local/add_run_info_to_multiqc_header/main'
include { SUBSAMPLE_FASTQ                } from '../modules/local/subsample/fastq/seqkit_head/main'
include { DOWNSAMPLED_BWA_METRICS        } from '../subworkflows/local/downsampled_bwa_metrics/main'
include { paramsSummaryMap               } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { removeKey ; partitioningActive ; appendExtension ; replaceExtension ; resolveRefPathsFromIndexDir ; resolveChrmPathsFromIndexDir
        ; resolveGiraffeIndexComponents ; makeFileChannel ; makeDirChannel ; makeFilePairChannel ; makeFastaChannel ; resolveAdditionalAlignmentMetricsRegions
        ; resolveVcReferences ; pairTumorNormalBams } from '../subworkflows/local/utils_xoos_workflow/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XOOS {
    take:
    ch_run_sheet // channel: run sheet read in from --input

    main:

    // Parameter validation

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Channel used to gather metric report JSON files from different steps across different
    // samples and runs to be used as input for the final metrics report generation step
    // during which we create one combined report for each run.
    // channel: [ val(run_meta), path(report) ]
    ch_run_metrics_report_json = channel.empty()
    ch_svc_metrics = channel.empty()

    // Custom MultiQC files that are generated by XOOS tools
    // channel: [ val(run_meta), path(multiqc_report_dir) ]
    ch_custom_multiqc_files_with_meta = channel.empty()

    // MultiQC files that are generated by modules supported by MultiQC (e.g. Picard)
    // channel: [ val(run_meta), path(multiqc_report_dir) ]
    ch_multiqc_files_with_meta = channel.empty()

    // Alignment metrics files for publishing
    // channel: [ val(meta), path(read_metrics), path(coverage), path(accuracy) ]
    ch_alignment_metrics_files = channel.empty()

    // GC bias metrics files for publishing
    // channel: [ val(meta), path(gc_bias) ]
    ch_gc_bias_metrics_files = channel.empty()

    // Hybrid selection (HS) metrics files for publishing
    // channel: [ val(meta), path(*_metrics) ]
    ch_hs_metrics_files = channel.empty()

    // Read collapser metrics files for publishing (summary_stats.tsv,
    // cluster_size_distribution_summary.tsv, cluster_size_distributions.tsv)
    // channel: [ val(meta), path(read_collapser_metrics_dir) ]
    ch_read_collapser_metrics_files = channel.empty()

    // Somalier outputs for publishing
    // channel: [ val(meta), path(pairs_tsv), path(samples_tsv), path(html) ]
    ch_somalier_output = channel.empty()

    // STEP 1: Route input runs by file_type into per-type channels
    //
    // The reference channel is built here (rather than in STEP 3) because
    // ROUTE_INPUTS converts pre-aligned CRAM inputs to BAM, which needs the
    // reference for decoding. Each CRAM is first checked for M5 tags by
    // VERIFY_CRAM so a CRAM lacking reference checksums fails fast.
    ch_reference = makeFastaChannel('reference', params.reference_fasta, params.reference_fai)

    ROUTE_INPUTS(ch_run_sheet, ch_reference)
    ch_basecall_files       = ROUTE_INPUTS.out.basecall_files
    ch_demultiplexed_fastqs = ROUTE_INPUTS.out.demultiplexed_fastqs
    ch_bams                 = ROUTE_INPUTS.out.bams
    // Keep the pre-aligned BAM/CRAM stream separate so mitochondrial BAM
    // post-processing can attach a filtered Giraffe graph to these samples.
    // ROUTE_INPUTS samples do not pass through ALIGN and therefore do not
    // produce a per-sample graph tuple.
    ch_pre_aligned_bams     = ch_bams
    ch_run_sample_names     = ROUTE_INPUTS.out.run_sample_names

    // STEP 2: Demultiplex basecall files

    // Note: Use value channel for demux resources since they should be the same for all runs
    ch_strand_detect_index = params.enable_demux_strand_detect.toBoolean()
        ? channel.value(file(params.demux_strand_detect_index, checkIfExists: true))
        : channel.value([])

    demux_out = DEMULTIPLEX(
        ch_basecall_files,
        ch_strand_detect_index,
    )
    ch_run_metrics_report_json = ch_run_metrics_report_json.concat(demux_out.report)
    ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.concat(demux_out.multiqc_files.map { meta, multiqc_report_dir -> [[id: meta.run_name, run_type: meta.run_type], multiqc_report_dir] })

    // Publish the raw per-run demux metrics directory (run_stats.tsv,
    // sample_stats.tsv, and read-length distribution tables) so consumers
    // have the unaggregated demux output alongside the aggregated report.
    ch_demux_metrics_files = demux_out.metrics_files

    // Combine pre-demultiplexed and newly demultiplexed FASTQ samples.
    // demultiplexed_fastqs contains only full reads when exclude_partial_reads
    // is enabled (all-split mode), so downstream steps never see partial reads.
    ch_all_fastqs = ch_demultiplexed_fastqs.concat(demux_out.demultiplexed_fastqs)

    // For publishing, use the complete sample directories (both full and
    // partial subdirectories) so that consumers of the published output
    // have access to all reads.
    ch_all_fastqs_for_publish = ch_demultiplexed_fastqs.concat(demux_out.demultiplexed_fastqs_all)

    // Optionally publish demultiplexed FASTQs
    ch_demuxed_fastqs_output = params.enable_publish_demuxed_fastqs.toBoolean()
        ? ch_all_fastqs_for_publish
        : channel.empty()

    // STEP 3: Alignment

    // Handle Subsampling of Demultiplexed FASTQs before alignment
    ch_pre_align_fastqs = ch_all_fastqs

    if (params.enable_subsample_demuxed_fastq.toBoolean()) {
        ch_fastq_to_subsample = ch_pre_align_fastqs.map { meta, fastqs ->
            def read_count = params.subsample_demuxed_fastq_read_count
            [ meta, read_count, fastqs ]
        }

        SUBSAMPLE_FASTQ(ch_fastq_to_subsample)
        ch_pre_align_fastqs = SUBSAMPLE_FASTQ.out.fastq.map { meta, fq -> [meta, [fq]] }
    }

    // The reference channel (ch_reference) is built in STEP 1 because
    // ROUTE_INPUTS needs it for CRAM->BAM conversion.

    // The post-read-collapser re-alignment in the consensus dedup path may use a
    // different aligner than the initial alignment. When dedup_aligner is set and
    // the consensus strategy is active it drives that re-alignment and produces
    // the final BAM; otherwise the single 'aligner' is used throughout.
    // Downstream steps that consume the final BAM (e.g. variant-caller model
    // selection) use final_aligner; steps operating on freshly aligned pre-dedup
    // BAMs use params.aligner.
    def final_aligner = (params.dedup_strategy == 'consensus' && params.dedup_aligner) ? params.dedup_aligner : params.aligner

    // Guard the pre-/post-read-collapser aligner combination. Two constraints:
    //   1. Giraffe aligners map against the linearized pangenome reference while
    //      the other aligners map against the standard reference, so mixing a
    //      Giraffe aligner with a non-Giraffe aligner across the read-collapser
    //      boundary would produce BAMs in different reference spaces.
    //   2. The two Giraffe aligners cannot be swapped for each other either:
    //      ch_giraffe_ref_paths (built below) is derived from params.aligner and
    //      pb_giraffe/vg_giraffe have different --ref-paths requirements
    //      (pb_giraffe strictly requires the .ref_paths_dict, vg_giraffe treats
    //      it as optional), so a Giraffe re-alignment must use the same Giraffe
    //      aligner as the initial alignment.
    // This is enforced at runtime rather than in the JSON schema because a
    // cross-field schema constraint on 'aligner' trips nf-core's default-parameter
    // schema linting (the default aligner is a Giraffe aligner).
    if (params.dedup_strategy == 'consensus' && params.dedup_aligner) {
        def giraffe_aligners = ['pb_giraffe', 'vg_giraffe']
        def aligner_is_giraffe = params.aligner in giraffe_aligners
        def dedup_aligner_is_giraffe = params.dedup_aligner in giraffe_aligners
        if (aligner_is_giraffe != dedup_aligner_is_giraffe) {
            error("Cannot mix Giraffe and non-Giraffe aligners across the read-collapser boundary: aligner='${params.aligner}', dedup_aligner='${params.dedup_aligner}'. Giraffe aligners (pb_giraffe, vg_giraffe) map against the linearized pangenome reference while the other aligners (bwa, minimap2, pb_fq2bam, pb_minimap2) map against the standard reference. Set both aligner and dedup_aligner to the same reference space (both Giraffe or both non-Giraffe).")
        }
        if (aligner_is_giraffe && params.aligner != params.dedup_aligner) {
            error("Cannot use different Giraffe aligners across the read-collapser boundary: aligner='${params.aligner}', dedup_aligner='${params.dedup_aligner}'. pb_giraffe and vg_giraffe have different --ref-paths requirements, and the reference paths used for the re-alignment are derived from aligner. Set dedup_aligner to the same Giraffe aligner as aligner.")
        }
    }

    // By default, recover chrM reads whenever mitochondrial calling is enabled
    // on a Giraffe-aligned BAM. An explicit false remains an opt-out.
    def enable_giraffe_mitochondrial_bam_post_process = params.enable_giraffe_mitochondrial_bam_post_process != null
        ? params.enable_giraffe_mitochondrial_bam_post_process.toBoolean()
        : (params.enable_mitochondrial_variant_caller.toBoolean() && params.aligner in ['pb_giraffe', 'vg_giraffe'])

    if (enable_giraffe_mitochondrial_bam_post_process) {
        if (!params.enable_mitochondrial_variant_caller.toBoolean() || params.variant_calling_mode != 'germline') {
            error("enable_giraffe_mitochondrial_bam_post_process requires germline mitochondrial variant calling. Set variant_calling_mode='germline' and enable_mitochondrial_variant_caller=true, or disable enable_giraffe_mitochondrial_bam_post_process.")
        }
        if (!(params.aligner in ['pb_giraffe', 'vg_giraffe'])) {
            error("enable_giraffe_mitochondrial_bam_post_process requires a Giraffe-aligned, NR-tagged BAM. Set aligner to 'pb_giraffe' or 'vg_giraffe', or disable enable_giraffe_mitochondrial_bam_post_process.")
        }
    }

    // Build the BWA index channel when either the initial aligner or the
    // post-read-collapser aligner needs it, so the consensus re-alignment can run
    // a BWA-based aligner even if the initial alignment did not.
    ch_bwa_idx = (params.aligner in ['pb_fq2bam', 'bwa'] || final_aligner in ['pb_fq2bam', 'bwa'])
        ? channel.value([[id: 'bwa_idx'], file(params.bwa_index_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])

    ch_giraffe_idx = params.aligner in ['pb_giraffe', 'vg_giraffe']
        ? channel.value(
            [
                [id: 'giraffe_idx'],
                params.enable_pangenome_personalization.toBoolean()
                    ? file(params.unfiltered_giraffe_index_dir, type: 'dir', checkIfExists: true)
                    : file(params.giraffe_index_dir, type: 'dir', checkIfExists: true),
            ]
        )
        : channel.value([[], [], []])

    ch_giraffe_chrm_paths = channel.value([])
    if (enable_giraffe_mitochondrial_bam_post_process
        && params.aligner in ['pb_giraffe', 'vg_giraffe']) {
        def chrm_paths_index_dir = params.enable_pangenome_personalization.toBoolean()
            ? params.unfiltered_giraffe_index_dir
            : params.giraffe_index_dir
        ch_giraffe_chrm_paths = channel.value(resolveChrmPathsFromIndexDir(chrm_paths_index_dir))
    }

    // --ref-paths for pbrun giraffe so the BAM @SQ header order matches the
    // reference used by downstream variant callers.
    // Without personalization: use the pre-built .dict derived from the linearized reference.
    // With personalization: use the pre-computed .paths file from the
    // unfiltered giraffe index directory (same prefix as the .gbz).
    if (params.aligner in ['pb_giraffe', 'vg_giraffe'] && params.enable_pangenome_personalization.toBoolean()) {
        ch_giraffe_ref_paths = channel.value(resolveRefPathsFromIndexDir(params.unfiltered_giraffe_index_dir))
    } else if (params.aligner == 'pb_giraffe' && params.giraffe_reference_linearized) {
        ch_giraffe_ref_paths = channel.value(file(replaceExtension(params.giraffe_reference_linearized, '.ref_paths_dict'), checkIfExists: true))
    } else if (params.aligner == 'vg_giraffe' && params.giraffe_reference_linearized) {
        // vg giraffe --ref-paths requires path names matching the GBZ graph; use the
        // sibling .ref_paths_dict only when it exists alongside the linearized reference.
        def vg_ref_paths_file = file(replaceExtension(params.giraffe_reference_linearized, '.ref_paths_dict'))
        ch_giraffe_ref_paths = vg_ref_paths_file.exists() ? channel.value(vg_ref_paths_file) : channel.value([])
    } else {
        ch_giraffe_ref_paths = channel.value([])
    }

    ch_alignment_input = ch_pre_align_fastqs

    if (params.enable_pangenome_personalization.toBoolean()) {
        // Check if sample-level metadata contains 'TN' and `group` keys to identify tumor-normal pairs
        ch_paired_samples = ch_pre_align_fastqs
            .filter { meta, _fastqs -> meta.containsKey('TN') && meta.containsKey('group') }
            .map { meta, fastqs -> [[group: meta.group], [meta, fastqs]] }
            .groupTuple(size: 2, remainder: true)
        // Concat any unpaired samples with the channel of samples that do not have the TN and group field
        // ch_unpaired_samples: [meta, fastqs]
        ch_unpaired_samples = ch_paired_samples
            .filter { _group, samples -> samples.size() == 1 }
            .map { _group, samples -> samples[0] }
            .mix(ch_pre_align_fastqs.filter { meta, _fastqs -> !(meta.containsKey('TN') && meta.containsKey('group')) })
        // ch_alignment_input: [meta, fastqs] where meta contains the field 'personalized_on'
        // which indicates which sample the personalized reference should be generated for
        ch_alignment_input = ch_paired_samples
            .filter { _group, samples -> samples.size() == 2 }
            .flatMap { group, samples ->
                def tumor = samples.find { meta, _fastqs -> meta.TN && meta.TN.toString().toUpperCase() == 'T' }
                def normal = samples.find { meta, _fastqs -> meta.TN && meta.TN.toString().toUpperCase() == 'N' }
                if (tumor && normal) {
                    // For the normal sample, set personalized_on to itself
                    normal[0].personalized_on = normal[1]
                    // For the tumor sample, set personalized_on to the normal sample
                    tumor[0].personalized_on = normal[1]
                    return [normal, tumor]
                }
                else {
                    log.warn("Could not identify tumor and normal FASTQ files for group ${group.group}. These samples will be treated as unpaired samples and the personalized reference will be generated using the sample's own FASTQ files.")
                    return samples.collect { meta, fastqs ->
                        [meta + [personalized_on: fastqs], fastqs]
                    }
                }
            }
            .mix(
                ch_unpaired_samples.map { meta, fastqs ->
                    // For unpaired samples, set personalized_on to itself
                    [meta + [personalized_on: fastqs], fastqs]
                }
            )
    }

    // Stamp meta.sample_id before alignment so every aligner can use it
    // for read group construction. Partitioning rewrites meta.id to
    // "SampleName_N", but sample_id stays as the true sample name.
    ch_alignment_input = ch_alignment_input.map { meta, reads ->
        [meta + [sample_id: meta.id], reads]
    }

    // Align all FASTQ samples and combine with pre-aligned BAMs
    ch_giraffe_graphs = channel.empty()
    if (partitioningActive()) {
        align_out = PARTITION_ALIGN_MERGE(
            params.aligner,
            params.alignment_partition_count.toInteger(),
            ch_alignment_input,
            ch_reference,
            ch_bwa_idx,
            ch_giraffe_idx,
            ch_giraffe_ref_paths,
        )
    } else {
        align_out = ALIGN(
            params.aligner,
            ch_alignment_input,
            ch_reference,
            ch_bwa_idx,
            ch_giraffe_idx,
            ch_giraffe_ref_paths,
        )
    }
    // Keep the Giraffe graph used for alignment
    // If a BAM is provided to the pipeline already aligned, the original Giraffe is unknown
    // and we fall back to using the filtered/base Giraffe graph.
    if (params.aligner in ['pb_giraffe', 'vg_giraffe']) {
        ch_giraffe_graphs = align_out.giraffe_graphs
            .map { meta, graph, zipcodes, minimizers, dist ->
                [removeKey(meta, 'personalized_on'), graph, zipcodes, minimizers, dist]
            }

        if (enable_giraffe_mitochondrial_bam_post_process) {
            // The graph used for a pre-aligned BAM is unknown. Use the
            // filtered/base Giraffe.
            ch_pre_aligned_giraffe_graphs = ch_pre_aligned_bams.map { meta, _bam, _bai ->
                def components = resolveGiraffeIndexComponents(params.giraffe_index_dir)
                [removeKey(meta, 'personalized_on'), components.gbz, components.zipcodes, components.minimizer, components.dist]
            }
            ch_giraffe_graphs = ch_giraffe_graphs.mix(ch_pre_aligned_giraffe_graphs)
        }
    }
    ch_bams = ch_bams.mix(
        align_out.bams.map { meta, bam, bai -> [removeKey(meta, 'personalized_on') + [bam_dedup_type: 'none', bam_secondary_rescued: false], bam, bai] }
    )

    // STEP 4: Build shared resource channels for metrics and variant calling

    ch_autosomes_bed                = makeFileChannel('autosomes_bed', params.autosomes_bed)
    ch_target_bed                   = makeFileChannel('target_bed', params.target_bed)
    ch_capture_bed                  = params.capture_bed
        ? channel.value([[id: 'capture_bed'], file(params.capture_bed, checkIfExists: true)])
        : channel.value([[id: 'capture_bed'], []])
    ch_high_confidence_bed_dir      = params.high_confidence_bed_dir
        ? channel.value([[id: 'high_confidence_bed_dir'], file(params.high_confidence_bed_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], [], [], []])
    ch_fallback_high_confidence_bed = makeFileChannel('fallback_high_confidence_bed', params.fallback_high_confidence_bed)

    def additional_bed_files = resolveAdditionalAlignmentMetricsRegions(
        params.alignment_metrics_additional_regions ?: [],
        params.alignment_metrics_additional_region_beds ?: [:],
    )

    ch_tfe_included_regions         = makeFileChannel('tfe_included_regions', params.tfe_included_regions_bed)
    ch_tfe_panel_of_normals         = makeFilePairChannel('tfe_panel_of_normals', params.tfe_panel_of_normals_tsv, params.tfe_panel_of_normal_tbi)
    ch_pop_af_vcf_with_tbi = params.pop_af_vcf
        ? channel.value([[id: 'pop_af_vcf'], file(params.pop_af_vcf, checkIfExists: true), file(params.pop_af_tbi, checkIfExists: true)])
        : channel.value([[], [], []])

    ch_contamination_vcf            = makeFileChannel('contamination_vcf', params.contamination_vcf)
    ch_contamination_included_regions = makeFileChannel('contamination_included_regions', params.contamination_included_regions_bed)
    ch_somalier_sites_vcf           = makeFileChannel('somalier_sites_vcf', params.somalier_sites_vcf)

    // Metrics reference: Giraffe has its own linearized reference created from extracting (linearizing) the
    // primary path of the graph reference. Use the linearized reference for metrics if aligned with
    // Giraffe (vg_giraffe or pb_giraffe); otherwise use the regular linear reference.
    def candidate_linearized_ref = params.enable_pangenome_personalization.toBoolean()
        ? params.unfiltered_giraffe_reference_linearized
        : params.giraffe_reference_linearized
    def reference_for_metrics = params.aligner in ['pb_giraffe', 'vg_giraffe'] && candidate_linearized_ref
        ? candidate_linearized_ref
        : params.reference_fasta
    ch_reference_for_metrics = reference_for_metrics
        ? channel.value([[id: 'reference_for_metrics'], file(reference_for_metrics, checkIfExists: true), file(appendExtension(reference_for_metrics, '.fai'), checkIfExists: true)])
        : channel.value([[], [], []])
    // Whether this run needs secondary-alignment rescue.
    def enable_secondary_rescue = (params.aligner in ['vg_giraffe', 'pb_giraffe'] && params.enable_copy_number_caller.toBoolean())

    // STEP 4.5: Downsampled BWA Metrics (optional, off by default)
    // Branches from ch_all_fastqs (full-depth, pre-subsample) to avoid
    // double-downsampling when enable_subsample_demuxed_fastq is also enabled.
    // Runs in parallel with the primary alignment path. Skipped for BAM inputs
    // since no FASTQs are available.
    if (params.enable_downsampled_bwa_metrics.toBoolean()) {
        ch_bwa_idx_for_metrics = channel.value([[id: 'bwa_idx'], file(params.bwa_index_dir, type: 'dir', checkIfExists: true)])
        // For downsampled BWA metrics, we align a subset of the reads using BWA (or fq2bam)
        // so we should always use the regular linear reference for metrics even if the main alignment was done with Giraffe.
        ch_reference_for_bwa_metrics = makeFastaChannel('reference_for_metrics', params.reference_fasta, params.reference_fai)

        DOWNSAMPLED_BWA_METRICS(
            ch_all_fastqs,
            ch_reference,
            ch_bwa_idx_for_metrics,
            ch_autosomes_bed,
            ch_target_bed,
            ch_capture_bed,
            ch_high_confidence_bed_dir,
            ch_fallback_high_confidence_bed,
            ch_reference_for_bwa_metrics,
            params.enable_alignment_metrics_autosomes.toBoolean(),
            params.enable_alignment_metrics_targets.toBoolean(),
            params.enable_alignment_metrics_high_confidence.toBoolean(),
            params.enable_gc_bias_metrics.toBoolean(),
            params.enable_hs_metrics.toBoolean(),
            params.downsampled_bwa_metrics_read_count,
            params.downsampled_bwa_metrics_aligner,
            additional_bed_files,
        )
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(DOWNSAMPLED_BWA_METRICS.out.report)
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
            DOWNSAMPLED_BWA_METRICS.out.custom_multiqc_files
                .map { meta, multiqc_dir -> [[id: meta.run_name, run_type: meta.run_type], multiqc_dir] }
        )
        ch_multiqc_files_with_meta = ch_multiqc_files_with_meta.mix(
            DOWNSAMPLED_BWA_METRICS.out.multiqc_files
                .map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
        )
        ch_alignment_metrics_files = ch_alignment_metrics_files.mix(DOWNSAMPLED_BWA_METRICS.out.alignment_metrics_files)
        ch_gc_bias_metrics_files = ch_gc_bias_metrics_files.mix(DOWNSAMPLED_BWA_METRICS.out.gc_bias_metrics_files)
        ch_hs_metrics_files = ch_hs_metrics_files.mix(DOWNSAMPLED_BWA_METRICS.out.hs_metrics_files)
    }

    // Compatible externally deduplicated BAMs retain their state through DEDUP.
    // Incompatible bam_type/dedup_strategy combinations fail during input routing.
    def enable_dedup = params.dedup_strategy != 'none'

    ch_raw_bams = ch_bams.filter { meta, _bam, _bai ->
        !(meta.bam_dedup_type in ['xoos_markdup', 'xoos_consensus'])
    }

    // STEP 5: Pre-dedup metrics (optional)
    // Calculate pre-dedup metrics if
    // - dedup is disabled and either subsampling is disabled (the pre-dedup BAM is the final BAM used for downstream analysis) or
    // - dedup is disabled and subsampling is enabled; pre-dedup metrics are explicitly enabled or
    // - dedup is enabled and pre-dedup metrics are explicitly enabled
    // If both dedup and subsampling are disabled, the pre-dedup BAM is the final BAM used for downstream analysis
    def calculate_metrics_on_original_bam =
        params.enable_pre_dedup_metrics.toBoolean() ||
        (!enable_dedup && (params.enable_pre_subsample_metrics.toBoolean() || !params.enable_subsample_bam.toBoolean()))
    if (calculate_metrics_on_original_bam) {
        ch_pre_dedup_bams = ch_raw_bams.filter { meta, _bam, _bai ->
            // Skip BAMs that have already been deduplicated externally.
            if (meta.bam_dedup_type != 'none') {
                log.warn("Skipping pre-dedup metrics for sample '${meta.id}': bam_type '${meta.bam_dedup_type}' is not eligible for pre-dedup metrics.")
                false
            }
            else {
                true
            }
        }
        PRE_DEDUP_METRICS(
            'original',
            ch_pre_dedup_bams,
            ch_autosomes_bed,
            ch_target_bed,
            ch_capture_bed,
            ch_high_confidence_bed_dir,
            ch_fallback_high_confidence_bed,
            ch_reference_for_metrics,
            params.enable_alignment_metrics_autosomes.toBoolean(),
            params.enable_alignment_metrics_targets.toBoolean(),
            params.enable_alignment_metrics_high_confidence.toBoolean(),
            params.enable_gc_bias_metrics.toBoolean(),
            params.enable_hs_metrics.toBoolean(),
            params.aligner,
            additional_bed_files,
        )
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(PRE_DEDUP_METRICS.out.report)
    }

    // STEP 6: Dedup (markdup or consensus)

    ch_pre_dedup_bam = channel.empty()
    if (enable_dedup) {
        if (params.enable_publish_pre_dedup_bam.toBoolean()) {
            ch_pre_dedup_bam = ch_bams
        }
        // Pass every BAM through DEDUP. Inputs that have already been deduped
        // are passed through unchanged.
        DEDUP(
            params.dedup_strategy,
            final_aligner,
            enable_secondary_rescue,
            ch_bams,
            ch_target_bed,
            ch_reference,
            ch_bwa_idx,
            ch_giraffe_idx,
            ch_giraffe_ref_paths,
        )
        ch_bams = DEDUP.out.bam
        if (params.dedup_strategy == 'consensus' && final_aligner in ['pb_giraffe', 'vg_giraffe']) {
            ch_giraffe_graphs = DEDUP.out.giraffe_graphs
                .map { meta, graph, zipcodes, minimizers, dist ->
                    [removeKey(meta, 'personalized_on'), graph, zipcodes, minimizers, dist]
                }
        }
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(DEDUP.out.report)

        // Publish the raw read_collapser metric TSVs (summary_stats.tsv,
        // cluster_size_distribution_summary.tsv, cluster_size_distributions.tsv).
        ch_read_collapser_metrics_files = ch_read_collapser_metrics_files.mix(DEDUP.out.read_collapser_metrics_dir)

        // Convert read_collapser metrics to MultiQC format
        XOOS_READ_COLLAPSER_METRICS_TO_MULTIQC(DEDUP.out.read_collapser_metrics_dir)
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
            XOOS_READ_COLLAPSER_METRICS_TO_MULTIQC.out.multiqc_report_dir
                .map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
        )
        // Calculate pre-subsample dedup metrics if
        // - subsampling is disabled (the deduped BAM is the final BAM used for downstream analysis)
        // - subsampling is enabled but pre-subsample metrics are explicitly enabled.
        if (params.enable_pre_subsample_metrics.toBoolean() || !params.enable_subsample_bam.toBoolean()) {
            DEDUP_METRICS(
                'dedup',
                ch_bams,
                ch_autosomes_bed,
                ch_target_bed,
                ch_capture_bed,
                ch_high_confidence_bed_dir,
                ch_fallback_high_confidence_bed,
                ch_reference_for_metrics,
                params.enable_alignment_metrics_autosomes.toBoolean(),
                params.enable_alignment_metrics_targets.toBoolean(),
                params.enable_alignment_metrics_high_confidence.toBoolean(),
                params.enable_gc_bias_metrics.toBoolean(),
                params.enable_hs_metrics.toBoolean(),
                params.aligner,
                additional_bed_files,
            )
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(DEDUP_METRICS.out.report)
        }
    }

    // STEP 7: Subsample BAM by coverage (optional)

    // Pre-subsample BAM (the BAM fed into subsampling). Only populated when
    // both subsampling and pre-subsample BAM publishing are enabled.
    ch_pre_subsample_bam = channel.empty()

    if (params.enable_subsample_bam.toBoolean()) {
        // Capture the pre-subsample BAM before ch_bams is overwritten by the
        // subsampled output, so it can optionally be published.
        if (params.enable_publish_pre_subsample_bam.toBoolean()) {
            ch_pre_subsample_bam = ch_bams
        }
        ch_bams = SUBSAMPLE_BAM(
            ch_bams,
            ch_reference,
            ch_autosomes_bed,
            params.subsample_bam_coverage_target,
            tuple('median', params.subsample_bam_metric_name),
        ).bam_subsampled
        // If subsampling is enabled, always calculate selected metrics on the subsampled BAM
        // as this is the final output BAM that will be used for downstream analysis.
        SUBSAMPLE_METRICS(
            enable_dedup ? 'dedup_subsample' : 'subsample',
            ch_bams,
            ch_autosomes_bed,
            ch_target_bed,
            ch_capture_bed,
            ch_high_confidence_bed_dir,
            ch_fallback_high_confidence_bed,
            ch_reference_for_metrics,
            params.enable_alignment_metrics_autosomes.toBoolean(),
            params.enable_alignment_metrics_targets.toBoolean(),
            params.enable_alignment_metrics_high_confidence.toBoolean(),
            params.enable_gc_bias_metrics.toBoolean(),
            params.enable_hs_metrics.toBoolean(),
            params.aligner,
            additional_bed_files,
        )
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(SUBSAMPLE_METRICS.out.report)
    }

    // STEP 7.5: Contamination Detection

    ch_contamination_output = channel.empty()

    if (params.enable_contamination_estimation.toBoolean()) {
        CONTAMINATION_DETECTION(
            ch_bams,
            ch_reference,
            ch_contamination_vcf,
            ch_contamination_included_regions,
        )
        ch_contamination_output = CONTAMINATION_DETECTION.out.results
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
            XOOS_PARSE_CONTAMINATION_METRICS(CONTAMINATION_DETECTION.out.results).report
        )
    }

    // STEP 7.75: Tumor Fraction Estimation

    ch_tfe_output = channel.empty()

    if (params.enable_tfe.toBoolean()) {
        TUMOR_FRACTION_ESTIMATION(
            ch_bams,
            ch_reference,
            ch_tfe_included_regions,
            ch_tfe_panel_of_normals
        )
        ch_tfe_output = TUMOR_FRACTION_ESTIMATION.out.results
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
            XOOS_PARSE_TUMOR_FRACTION_METRICS(TUMOR_FRACTION_ESTIMATION.out.results).report
        )
    }

    // STEP 8: Variant calling

    ch_vcf = channel.empty()
    ch_raw_vc_vcf = channel.empty()
    ch_mutect_stats = channel.empty()
    ch_svc_pretrained_vcf = channel.empty()
    ch_svc_retrained_vcf = channel.empty()
    ch_svc_mitochondrial_vcf = channel.empty()
    ch_vcfeval_output = channel.empty()
    ch_cnc_output = channel.empty()
    ch_msi_output = channel.empty()
    ch_germline_filtering_output = channel.empty()
    ch_str_output = channel.empty()

    // Build variant calling reference channels (linearized ref / subset ref when using Giraffe)
    def vcRefs = resolveVcReferences(params, ch_reference)
    ch_vc_reference                  = vcRefs.reference
    ch_vc_reference_dict             = vcRefs.referenceDict
    ch_splitintervals_reference      = vcRefs.splitintervalsReference
    ch_splitintervals_reference_dict = vcRefs.splitintervalsReferenceDict

    ch_germline_svc_truth_vcf_dir   = makeDirChannel('truth_vcf_dir', params.germline_small_variant_truth_vcf_dir)
    ch_germline_svc_train_bed_dir   = makeDirChannel('train_bed_dir', params.germline_small_variant_train_bed_dir)
    ch_germline_svc_test_bed_dir    = makeDirChannel('test_bed_dir', params.germline_small_variant_test_bed_dir)
    ch_germline_svc_eval_bed_dir    = makeDirChannel('eval_bed_dir', params.germline_small_variant_eval_bed_dir)
    ch_somatic_tn_svc_truth_vcf_dir = makeDirChannel('somatic_tn_truth_vcf_dir', params.somatic_tumor_normal_small_variant_truth_vcf_dir)
    ch_somatic_tn_svc_eval_bed_dir  = makeDirChannel('somatic_tn_eval_bed_dir', params.somatic_tumor_normal_small_variant_eval_bed_dir)
    ch_somatic_to_svc_truth_vcf_dir = makeDirChannel('somatic_to_truth_vcf_dir', params.svc_somatic_to_truth_vcf_dir)
    ch_somatic_to_svc_eval_bed_dir  = makeDirChannel('somatic_to_eval_bed_dir', params.svc_somatic_to_eval_bed_dir)
    ch_svc_germline_models          = makeFilePairChannel('svc_models_germline', params.svc_model_germline_snv, params.svc_model_germline_indel)
    ch_svc_somatic_tn_models        = makeFileChannel('svc_model_somatic_tumor_normal', params.svc_model_somatic_tumor_normal)

    if (params.svc_panel_of_normals_vcf && !params.svc_panel_of_normals_tbi) {
        error('svc_panel_of_normals_vcf requires svc_panel_of_normals_tbi. Provide both files, or omit both parameters.')
    }
    if (params.svc_panel_of_normals_tbi && !params.svc_panel_of_normals_vcf) {
        error('svc_panel_of_normals_tbi requires svc_panel_of_normals_vcf. Provide both files, or omit both parameters.')
    }

    ch_panel_of_normals_vcf_with_tbi = makeFilePairChannel('panel_of_normals', params.svc_panel_of_normals_vcf, params.svc_panel_of_normals_tbi)
    ch_svc_somatic_to_model         = makeFileChannel('svc_model_somatic_tumor_only', params.svc_model_somatic_tumor_only)
    ch_par_beds = params.chr_x_par_bed
        ? channel.value([[id: 'chr_x_par_bed'], file(params.chr_x_par_bed, checkIfExists: true), file(params.chr_y_par_bed, checkIfExists: true)])
        : channel.value([[], []])
    ch_mappability_bw               = makeFileChannel('mappability_bw', params.mappability_bw)
    ch_cnv_seed_file                = makeFileChannel('cnv_seed_file', params.cnv_seed_seg)
    ch_cnv_blocklist_bed            = makeFileChannel('cnv_blocklist_bed', params.cnv_blocklist_bed)
    ch_cnv_tumor_only_blocklist_bed = makeFileChannel('cnc_somatic_to_blocklist_bed', params.cnc_somatic_to_blocklist_bed)
    ch_cnv_snp_blocklist_bed        = makeFileChannel('cnc_somatic_to_snp_blocklist_bed', params.cnc_somatic_to_snp_blocklist_bed)
    ch_repeat_catalog               = makeFileChannel('repeat_catalog', params.repeat_catalog)
    ch_msi_sites_thresholds_tsv     = makeFileChannel('msi_sites_thresholds_tsv', params.msi_sites_thresholds_tsv)
    ch_msi_background_threshold_json = makeFileChannel('msi_background_threshold_json', params.msi_background_threshold_json)

    def reference_sdf = params.reference_sdf ?: null

    // Enable variant calling partitioning if the option is enabled and the partition count for the
    // chosen variant caller is greater than 1.
    def use_cpu_variant_callers =
        (params.variant_calling_mode in ['somatic_tn', 'somatic_to'] && params.use_gatk_mutect2.toBoolean()) ||
        (params.variant_calling_mode == 'germline' && params.use_gatk_haplotype_caller.toBoolean())
    def vc_partition_count = use_cpu_variant_callers
        ? params.cpu_variant_calling_partition_count.toInteger()
        : params.gpu_variant_calling_partition_count.toInteger()
    def vc_partitioning_enabled = vc_partition_count > 1 && params.enable_variant_calling_partitioning.toBoolean()

    if (params.variant_calling_mode == 'somatic_tn') {
        ch_tn_bam_pairs = pairTumorNormalBams(ch_bams)

        // Run somalier T/N pair verification if enabled
        if (params.enable_somatic_tn_pair_verification.toBoolean()) {
            SOMATIC_TN_PAIR_VERIFICATION(
                ch_tn_bam_pairs,
                ch_reference,
                ch_somalier_sites_vcf,
                params.somatic_tn_pair_mismatch_action,
                params.somatic_tn_pair_concordance_cutoff as Double
            )
            ch_tn_bam_pairs = SOMATIC_TN_PAIR_VERIFICATION.out.verified_pairs

            ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
                SOMATIC_TN_PAIR_VERIFICATION.out.multiqc_files
            )
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
                SOMATIC_TN_PAIR_VERIFICATION.out.report
            )

            ch_somalier_output = ch_somalier_output.mix(
                SOMATIC_TN_PAIR_VERIFICATION.out.somalier_output
            )
        }

        CALL_SOMATIC_VARIANTS_TN(
            ch_tn_bam_pairs,
            ch_vc_reference,
            ch_vc_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_target_bed,
            ch_pop_af_vcf_with_tbi,
            ch_svc_somatic_tn_models,
            ch_panel_of_normals_vcf_with_tbi,
            ch_somatic_tn_svc_truth_vcf_dir,
            ch_somatic_tn_svc_eval_bed_dir,
            ch_mappability_bw,
            ch_cnv_seed_file,
            ch_cnv_blocklist_bed,
            params.enable_small_variant_caller.toBoolean(),
            params.use_gatk_mutect2.toBoolean(),
            params.enable_mutect_stats.toBoolean(),
            params.enable_vcfeval.toBoolean(),
            params.enable_copy_number_caller.toBoolean(),
            params.enable_merge_vcf.toBoolean(),
            reference_sdf,
            vc_partitioning_enabled,
            vc_partition_count,
            final_aligner,
            ch_msi_sites_thresholds_tsv,
            params.enable_msi_caller.toBoolean(),
        )
        // Only publish through the generic vcf target when the callers were
        // actually merged into a single VCF. For somatic_tn merging is always
        // disabled (enforced by schema validation), so merged_vcf carries the
        // raw per-caller VCFs, which are already published by their dedicated
        // targets (svc_pretrained_vcf, cnc_output). Feeding them to the vcf
        // target as well would duplicate the small variant caller VCF and
        // scatter the copy number caller index as an orphaned .tbi.
        ch_vcf = params.enable_merge_vcf.toBoolean()
            ? CALL_SOMATIC_VARIANTS_TN.out.merged_vcf
            : channel.empty()
        ch_raw_vc_vcf = CALL_SOMATIC_VARIANTS_TN.out.raw_vc_vcf
            .map { meta, vcf, tbi -> [meta + [raw_variant_caller: 'mutect2'], vcf, tbi] }
        ch_mutect_stats = CALL_SOMATIC_VARIANTS_TN.out.mutect_stats
        ch_svc_pretrained_vcf = CALL_SOMATIC_VARIANTS_TN.out.svc_vcf
        ch_vcfeval_output = CALL_SOMATIC_VARIANTS_TN.out.small_variant_vcfeval_output
        ch_svc_metrics = ch_svc_metrics.mix(CALL_SOMATIC_VARIANTS_TN.out.svc_metrics)
        ch_cnc_output = CALL_SOMATIC_VARIANTS_TN.out.copy_number_caller_output
        ch_msi_output = CALL_SOMATIC_VARIANTS_TN.out.msi_caller_output
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(CALL_SOMATIC_VARIANTS_TN.out.small_variant_vcfeval_report)

        if (params.enable_copy_number_caller.toBoolean()) {
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
                XOOS_PARSE_CNC_METRICS_SOMATIC(CALL_SOMATIC_VARIANTS_TN.out.copy_number_caller_output).report
            )
        }
    }

    if (params.variant_calling_mode == 'somatic_to') {
        // Tumor-only: every sample is processed on its own, no TN pairing or
        // group column is required. Drop any pairing-related meta keys so the
        // downstream somatic logic treats each sample as a standalone tumor.
        ch_tumor_only_bams = ch_bams
            .map { meta, bam, bai ->
                def m = meta.clone()
                m.remove('TN')
                m.remove('group')
                m.remove('bam_dedup_type')
                m.remove('bam_secondary_rescued')
                m.sample_type = meta.sample_type ?: 'ffpe'
                // Return [meta, [tumor_bam], [tumor_bai]] to match the shape the
                // Mutect subworkflows expect (list of BAMs, list of indexes).
                return [m, [bam], [bai]]
            }

        CALL_SOMATIC_VARIANTS_TO(
            ch_tumor_only_bams,
            ch_vc_reference,
            ch_vc_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_target_bed,
            ch_pop_af_vcf_with_tbi,
            ch_svc_somatic_to_model,
            ch_panel_of_normals_vcf_with_tbi,
            ch_somatic_to_svc_truth_vcf_dir,
            ch_somatic_to_svc_eval_bed_dir,
            ch_mappability_bw,
            ch_cnv_seed_file,
            ch_cnv_tumor_only_blocklist_bed,
            ch_cnv_snp_blocklist_bed,
            params.enable_small_variant_caller.toBoolean(),
            params.use_gatk_mutect2.toBoolean(),
            params.enable_mutect_stats.toBoolean(),
            params.enable_vcfeval.toBoolean(),
            params.enable_copy_number_caller.toBoolean(),
            params.enable_germline_filtering.toBoolean(),
            params.enable_merge_vcf.toBoolean(),
            reference_sdf,
            vc_partitioning_enabled,
            vc_partition_count,
            ch_msi_background_threshold_json,
            params.enable_msi_caller.toBoolean(),
        )
        // Tumor-only merges by default: the small-variant contribution (the
        // germline-filtered VCF when enable_germline_filtering is set,
        // otherwise the small-variant VCF) is concatenated with the copy-number
        // VCF and, when enabled, the MSI VCF. The merged result is routed
        // through the generic vcf target. When merging is disabled, per-caller
        // VCFs are published by their dedicated targets.
        ch_vcf = params.enable_merge_vcf.toBoolean()
            ? CALL_SOMATIC_VARIANTS_TO.out.merged_vcf
            : channel.empty()
        ch_raw_vc_vcf = CALL_SOMATIC_VARIANTS_TO.out.raw_vc_vcf
            .map { meta, vcf, tbi -> [meta + [raw_variant_caller: 'mutect2'], vcf, tbi] }
        ch_mutect_stats = CALL_SOMATIC_VARIANTS_TO.out.mutect_stats
        // When germline filtering is enabled its VCF supersedes the raw SVC VCF for
        // publishing — this covers both the merge and no-merge paths.
        ch_svc_pretrained_vcf = params.enable_germline_filtering.toBoolean()
            ? CALL_SOMATIC_VARIANTS_TO.out.germline_filtered_vcf
            : CALL_SOMATIC_VARIANTS_TO.out.svc_vcf
        ch_vcfeval_output = CALL_SOMATIC_VARIANTS_TO.out.small_variant_vcfeval_output
        ch_svc_metrics = ch_svc_metrics.mix(CALL_SOMATIC_VARIANTS_TO.out.svc_metrics)
        ch_cnc_output = CALL_SOMATIC_VARIANTS_TO.out.copy_number_caller_output
        ch_germline_filtering_output = CALL_SOMATIC_VARIANTS_TO.out.germline_filtering_output
        ch_msi_output = CALL_SOMATIC_VARIANTS_TO.out.msi_caller_output
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(CALL_SOMATIC_VARIANTS_TO.out.small_variant_vcfeval_report)

        if (params.enable_copy_number_caller.toBoolean()) {
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
                XOOS_PARSE_CNC_METRICS_SOMATIC_TO(CALL_SOMATIC_VARIANTS_TO.out.copy_number_caller_output).report
            )
        }
    }

    if (params.variant_calling_mode == 'germline') {
        CALL_GERMLINE_VARIANTS(
            ch_bams,
            ch_vc_reference,
            ch_vc_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_target_bed,
            ch_par_beds,
            ch_pop_af_vcf_with_tbi,
            ch_svc_germline_models,
            ch_germline_svc_truth_vcf_dir,
            ch_germline_svc_train_bed_dir,
            ch_germline_svc_test_bed_dir,
            ch_germline_svc_eval_bed_dir,
            ch_mappability_bw,
            ch_cnv_seed_file,
            ch_cnv_blocklist_bed,
            ch_repeat_catalog,
            params.enable_small_variant_caller.toBoolean(),
            params.use_gatk_haplotype_caller.toBoolean(),
            final_aligner,
            params.enable_small_variant_caller_retrain.toBoolean(),
            params.enable_vcfeval.toBoolean(),
            params.enable_copy_number_caller.toBoolean(),
            params.enable_copy_number_caller_baf.toBoolean(),
            params.enable_str_caller.toBoolean(),
            params.enable_mitochondrial_variant_caller.toBoolean(),
            params.use_gatk_mutect2_mitochondrial.toBoolean(),
            params.enable_merge_vcf.toBoolean(),
            reference_sdf,
            vc_partitioning_enabled,
            vc_partition_count,
            ch_giraffe_graphs,
            ch_giraffe_chrm_paths,
        )
        ch_vcf = CALL_GERMLINE_VARIANTS.out.merged_vcf
        ch_raw_vc_vcf = CALL_GERMLINE_VARIANTS.out.raw_vc_vcf
            .map { meta, vcf, tbi -> [meta + [raw_variant_caller: 'haplotype_caller'], vcf, tbi] }
        ch_svc_pretrained_vcf = CALL_GERMLINE_VARIANTS.out.svc_pretrained_vcf
        ch_svc_retrained_vcf = CALL_GERMLINE_VARIANTS.out.svc_retrained_vcf
        ch_svc_mitochondrial_vcf = CALL_GERMLINE_VARIANTS.out.mitochondrial_vcf
        ch_vcfeval_output = CALL_GERMLINE_VARIANTS.out.small_variant_vcfeval_output
        ch_svc_metrics = ch_svc_metrics.mix(CALL_GERMLINE_VARIANTS.out.svc_metrics)
        ch_cnc_output = CALL_GERMLINE_VARIANTS.out.copy_number_caller_output
        ch_str_output = CALL_GERMLINE_VARIANTS.out.str_caller_output
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(CALL_GERMLINE_VARIANTS.out.small_variant_vcfeval_report)

        if (params.enable_str_caller.toBoolean()) {
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
                XOOS_PARSE_STR_METRICS(CALL_GERMLINE_VARIANTS.out.str_caller_output).report
            )
        }
        if (params.enable_copy_number_caller.toBoolean()) {
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
                XOOS_PARSE_CNC_METRICS_GERMLINE(CALL_GERMLINE_VARIANTS.out.copy_number_caller_output).report
            )
        }
    }

    // STEP 9: Metrics reporting

    // Convert SVC filter_variants TSV metrics to the per-sample JSON report format
    // consumed by the run-level report generator.
    ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
        XOOS_PARSE_SVC_METRICS(
            ch_svc_metrics.map { meta, metrics -> [meta, metrics, meta.svc_source] }
        ).report
    )

    // We only need sample name and run name for the metrics report
    // so we will extract these from the meta map to avoid any issues
    // with groupTuple. Use remainder:true so partial groups from ignored
    // task failures are emitted rather than causing a deadlock.
    def ch_grouped_metrics = ch_run_metrics_report_json
        .map { meta, report_path -> [[id: meta.run_name, run_type: meta.run_type], report_path] }
        .groupTuple(remainder: true)
        .map { meta, reports -> [meta.id, meta, reports] }
        .combine(ch_run_sample_names, by: 0)
        .map { _run_name, meta, reports, sample_names -> [meta + [sample_name: sample_names], reports] }

    ch_combined_metrics_report = XOOS_GENERATE_REPORT(ch_grouped_metrics).report

    ch_combined_metrics_report_summary = params.metrics_report_summary_patterns
        ? XOOS_GENERATE_REPORT_SUMMARY(ch_grouped_metrics, params.metrics_report_summary_patterns).summary_report
        : channel.empty()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        // remainder: true has no effect here because size is unspecified (Nextflow
        // GroupTupleOp only applies remainder to size-bounded groups). The topic
        // channel closing on session completion or abort is what allows this
        // groupTuple to emit. The root-cause fix for the abort-hang is ensuring
        // all topic publishers succeed (see VG_GIRAFFE and KMC eval fixes), so
        // the topic closes before the versions pipeline can block shutdown.
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'xoos_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    // The metrics reported should be from the final BAM used for downstream analysis
    if (params.enable_subsample_bam.toBoolean()) {
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(SUBSAMPLE_METRICS.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_multiqc_files_with_meta = ch_multiqc_files_with_meta.mix(SUBSAMPLE_METRICS.out.multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_alignment_metrics_files = ch_alignment_metrics_files.mix(SUBSAMPLE_METRICS.out.alignment_metrics_files)
        ch_gc_bias_metrics_files = ch_gc_bias_metrics_files.mix(SUBSAMPLE_METRICS.out.gc_bias_metrics_files)
        ch_hs_metrics_files = ch_hs_metrics_files.mix(SUBSAMPLE_METRICS.out.hs_metrics_files)
    }
    else if (enable_dedup) {
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(DEDUP_METRICS.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_multiqc_files_with_meta = ch_multiqc_files_with_meta.mix(DEDUP_METRICS.out.multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_alignment_metrics_files = ch_alignment_metrics_files.mix(DEDUP_METRICS.out.alignment_metrics_files)
        ch_gc_bias_metrics_files = ch_gc_bias_metrics_files.mix(DEDUP_METRICS.out.gc_bias_metrics_files)
        ch_hs_metrics_files = ch_hs_metrics_files.mix(DEDUP_METRICS.out.hs_metrics_files)
    }
    else if (!params.enable_subsample_bam.toBoolean()) {
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(PRE_DEDUP_METRICS.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_multiqc_files_with_meta = ch_multiqc_files_with_meta.mix(PRE_DEDUP_METRICS.out.multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_alignment_metrics_files = ch_alignment_metrics_files.mix(PRE_DEDUP_METRICS.out.alignment_metrics_files)
        ch_gc_bias_metrics_files = ch_gc_bias_metrics_files.mix(PRE_DEDUP_METRICS.out.gc_bias_metrics_files)
        ch_hs_metrics_files = ch_hs_metrics_files.mix(PRE_DEDUP_METRICS.out.hs_metrics_files)
    }

    if (params.enable_small_variant_caller.toBoolean() && params.enable_vcfeval.toBoolean()) {
        if (params.variant_calling_mode == 'germline') {
            ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(CALL_GERMLINE_VARIANTS.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        }
        else if (params.variant_calling_mode == 'somatic_tn') {
            ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(CALL_SOMATIC_VARIANTS_TN.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        }
    }

    // Contamination detection metrics
    if (params.enable_contamination_estimation.toBoolean()) {
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
            CONTAMINATION_DETECTION.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
        )
    }

    // STR caller metrics (germline only)
    if (params.variant_calling_mode == 'germline' && params.enable_str_caller.toBoolean()) {
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
            XOOS_STR_METRICS_TO_MULTIQC(CALL_GERMLINE_VARIANTS.out.str_caller_output).multiqc_report_dir
                .map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
        )
    }

    // Copy number caller metrics
    if (params.enable_copy_number_caller.toBoolean()) {
        if (params.variant_calling_mode == 'germline') {
            ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
                XOOS_CNC_METRICS_TO_MULTIQC_GERMLINE(CALL_GERMLINE_VARIANTS.out.copy_number_caller_output).multiqc_report_dir
                    .map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
            )
        }
        if (params.variant_calling_mode == 'somatic_tn') {
            ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
                XOOS_CNC_METRICS_TO_MULTIQC_SOMATIC(CALL_SOMATIC_VARIANTS_TN.out.copy_number_caller_output).multiqc_report_dir
                    .map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
            )
        }
        if (params.variant_calling_mode == 'somatic_to') {
            ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
                XOOS_CNC_METRICS_TO_MULTIQC_SOMATIC_TO(CALL_SOMATIC_VARIANTS_TO.out.copy_number_caller_output).multiqc_report_dir
                    .map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
            )
        }
    }

    // Tumor fraction estimation metrics
    if (params.enable_tfe.toBoolean()) {
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(
            TUMOR_FRACTION_ESTIMATION.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] }
        )
    }

    // channel: [ val(run_meta), [path(multiqc_report_dir), ...] ]
    // Use remainder:true so partial groups from ignored task failures are emitted
    // rather than causing a deadlock.
    //
    // NOTE: the meta here is deliberately kept as {id, run_type} (no sample_name
    // enrichment). It must match the meta shape of ch_multiqc_files_with_meta
    // (the Picard channel) so that the groupTuple below unifies both into a
    // single group per run. Enriching with sample_name here would create a
    // distinct meta map, splitting the group and spawning two MULTIQC tasks that
    // each see only half the inputs (one Picard-only, one XOOS-only), with the
    // last to publish overwriting the other. sample_name is added once, after
    // the mix, on the combined channel below.
    ch_multiqc_files = XOOS_COMBINE_MULTIQC_REPORTS(
        ch_custom_multiqc_files_with_meta.groupTuple(remainder: true)
    ).multiqc_report_dir

    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_files_with_meta)
    //
    // MODULE: MultiQC
    //
    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_global_multiqc_files = ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')


    ch_global_multiqc_files = ch_global_multiqc_files.mix(ch_collated_versions)

    // Enrich per-run MultiQC files with the sample_name list so nf-report can
    // attribute global task failures to all affected samples. Use remainder:true
    // so partial groups from ignored task failures are still emitted.
    //
    // NOTE: global MultiQC files (workflow summary, software versions, methods
    // description) are deliberately NOT mixed into ch_multiqc_files here — that
    // would create a cyclic dependency between ADD_RUN_INFO_TO_MULTIQC_HEADER's
    // input and its own `topic: versions` output (via ch_global_multiqc_files →
    // ch_collated_versions → channel.topic("versions")), causing a permanent
    // hang at the end of the run. Globals are joined into MULTIQC's input
    // separately below; MULTIQC itself does not emit to the versions topic
    // (see modules/nf-core/multiqc/main.nf comment) so it can safely depend on
    // ch_global_multiqc_files.
    ch_multiqc_files = ch_multiqc_files
        .groupTuple(remainder: true)
        .map { meta, dirs -> [meta.id, meta, dirs] }
        .combine(ch_run_sample_names, by: 0)
        .map { _run_name, meta, dirs, sample_names -> [meta + [sample_name: sample_names], dirs] }

    // Render run-level info templates. This will be displayed in the header of the MultiQC report for each run.
    // ADD_RUN_INFO emits to `topic: versions`, so its input must NOT transitively
    // depend on the versions topic — see the note above.
    ch_multiqc_config_with_header = ADD_RUN_INFO_TO_MULTIQC_HEADER(
        ch_multiqc_files.map { meta, _multiqc_report_dir -> [meta, file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)] }
    ).multiqc_config_with_run_info

    // Collect the global multiqc files into a value channel that emits the
    // list as a single tuple element. The extra `.map { [it] }` wrapping is
    // important: without it, `combine` flattens the List<Path> into individual
    // tuple elements, breaking the destructuring closure below.
    ch_global_multiqc_files_list = ch_global_multiqc_files.collect().map { [it] }

    MULTIQC(
        ch_multiqc_files
            .join(ch_multiqc_config_with_header)
            .combine(ch_global_multiqc_files_list)
            .map { run_meta, multiqc_report_dirs, multiqc_config_with_run_info, global_files ->
                [
                    run_meta,
                    (multiqc_report_dirs instanceof List ? multiqc_report_dirs : [multiqc_report_dirs]) + global_files,
                    multiqc_config_with_run_info,
                    params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
    )

    emit:
    multiqc_report          = MULTIQC.out.report // channel: [ val(run_meta), path(multiqc_report.html) ]
    bam                     = ch_bams // channel: [ [meta], path(bam), path(bai) ]
    pre_dedup_bam           = ch_pre_dedup_bam // channel: [ [meta], path(bam), path(bai) ]
    pre_subsample_bam       = ch_pre_subsample_bam // channel: [ [meta], path(bam), path(bai) ]
    vcf                     = ch_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    raw_vc_vcf              = ch_raw_vc_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    mutect_stats            = ch_mutect_stats // channel: [ [meta], path(stats) ]
    svc_pretrained_vcf      = ch_svc_pretrained_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    svc_retrained_vcf       = ch_svc_retrained_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    svc_mitochondrial_vcf   = ch_svc_mitochondrial_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    vcfeval_output          = ch_vcfeval_output // channel: [ [meta], path(vcfeval_dir) ]
    demuxed_fastqs          = ch_demuxed_fastqs_output // channel: [ [meta], path(sample_fastq_dir) ]
    demux_metrics_files     = ch_demux_metrics_files // channel: [ [meta], path(metrics_dir) ]
    cnc_output              = ch_cnc_output // channel: [ [meta], path(copy_number_caller_output) ]
    msi_output              = ch_msi_output // channel: [ [meta], path(msi_caller_output) ]
    germline_filtering_output = ch_germline_filtering_output // channel: [ [meta], path(prediction.csv) ]
    str_output              = ch_str_output // channel: [ [meta], path(str_caller_output) ]
    tfe_output              = ch_tfe_output // channel: [ [meta], path(tfe_results_dir) ]
    contamination_output    = ch_contamination_output // channel: [ [meta], path(contamination_results_dir) ]
    alignment_metrics_files = ch_alignment_metrics_files // channel: [ [meta], path(read_metrics), path(coverage), path(accuracy) ]
    gc_bias_metrics_files   = ch_gc_bias_metrics_files // channel: [ [meta], path(gc_bias) ]
    hs_metrics_files        = ch_hs_metrics_files // channel: [ [meta], path(*_metrics) ]
    read_collapser_metrics_files = ch_read_collapser_metrics_files // channel: [ [meta], path(read_collapser_metrics_dir) ]
    somalier_output         = ch_somalier_output // channel: [ [meta], path(pairs_tsv), path(samples_tsv), path(html) ]
    combined_report         = ch_combined_metrics_report // channel: [ val(run_meta), path(report) ]
    combined_report_summary = ch_combined_metrics_report_summary // channel: [ val(run_meta), path(report_summary) ]
}
