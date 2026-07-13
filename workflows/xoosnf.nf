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
include { SAMTOOLS_REHEADER              } from '../modules/nf-core/samtools/reheader/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_REHEADER } from '../modules/nf-core/samtools/index/main'
include { METRICS as PRE_DEDUP_METRICS ; METRICS as DEDUP_METRICS ; METRICS as SUBSAMPLE_METRICS } from '../subworkflows/local/metrics/main'
include { CONTAMINATION_DETECTION        } from '../subworkflows/local/metrics/contamination_detection.nf'
include { TUMOR_FRACTION_ESTIMATION      } from '../subworkflows/local/tumor_fraction_estimation/main'
include { CALL_GERMLINE_VARIANTS         } from '../subworkflows/local/call_germline_variants/main'
include { CALL_SOMATIC_VARIANTS_TN       } from '../subworkflows/local/call_somatic_variants_tn/main'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { XOOS_GENERATE_REPORT           } from '../modules/local/xoos_generate_metrics_report/generate_report/main'
include { XOOS_GENERATE_REPORT_SUMMARY   } from '../modules/local/xoos_generate_metrics_report/generate_report_summary/main'
include { XOOS_PARSE_CONTAMINATION_METRICS                                   } from '../modules/local/xoos_generate_metrics_report/parse_contamination_metrics/main'
include { XOOS_PARSE_TUMOR_FRACTION_METRICS                                   } from '../modules/local/xoos_generate_metrics_report/parse_tumor_fraction_metrics/main'
include { XOOS_PARSE_STR_METRICS                                               } from '../modules/local/xoos_generate_metrics_report/parse_str_metrics/main'
include { XOOS_PARSE_CNC_METRICS as XOOS_PARSE_CNC_METRICS_SOMATIC            } from '../modules/local/xoos_generate_metrics_report/parse_cnc_metrics/main'
include { XOOS_PARSE_CNC_METRICS as XOOS_PARSE_CNC_METRICS_GERMLINE           } from '../modules/local/xoos_generate_metrics_report/parse_cnc_metrics/main'
include { XOOS_COMBINE_MULTIQC_REPORTS   } from '../modules/local/xoos_metrics_to_multiqc/combine_reports/main'
include { XOOS_STR_METRICS_TO_MULTIQC                                        } from '../modules/local/xoos_metrics_to_multiqc/str_metrics/main'
include { XOOS_CNC_METRICS_TO_MULTIQC as XOOS_CNC_METRICS_TO_MULTIQC_SOMATIC  } from '../modules/local/xoos_metrics_to_multiqc/cnc_metrics/main'
include { XOOS_CNC_METRICS_TO_MULTIQC as XOOS_CNC_METRICS_TO_MULTIQC_GERMLINE } from '../modules/local/xoos_metrics_to_multiqc/cnc_metrics/main'
include { XOOS_READ_COLLAPSER_METRICS_TO_MULTIQC  } from '../modules/local/xoos_metrics_to_multiqc/read_collapser_metrics/main'
include { ADD_RUN_INFO_TO_MULTIQC_HEADER } from '../modules/local/add_run_info_to_multiqc_header/main'
include { SUBSAMPLE_BASECALL_RUNS        } from '../subworkflows/local/subsample_basecall_runs/main'
include { SUBSAMPLE_FASTQ                } from '../modules/local/subsample/fastq/main'
include { paramsSummaryMap ; samplesheetToList } from 'plugin/nf-schema'
include { paramsSummaryMultiqc           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { convertSamplesheet ; validateSamplesheet } from '../subworkflows/local/utils_samplesheet/main'

/*
 * Helper: remove a key from a map and return the new map
 */
def removeKey(map, key) {
    def new_map = map.clone()
    new_map.remove(key)
    return new_map
}

/*
 * Helper: whether alignment partitioning is effectively active.
 * Partitioning only takes effect when it is enabled AND the partition count is
 * greater than 1 (a count of 1 is a no-op single partition). Used to decide
 * whether per-partition duplicate marking was disabled (so the merged BAM still
 * needs marking) and how freshly aligned BAMs are stamped.
 */
def partitioningActive() {
    return params.enable_alignment_partitioning.toBoolean() && params.alignment_partition_count.toInteger() > 1
}

/*
 * Helper: read a samplesheet CSV and return list of [meta_with_id, row_map]
 */
def readSamplesheet(samplesheet_path) {
    def samplesheet_file = file(samplesheet_path, checkIfExists: true)
    // Use samplesheetToList to validate the samplesheet against the schema
    def _samplesheet = samplesheetToList(samplesheet_file, "${projectDir}/assets/schema_samplesheet.json")
    return samplesheet_file.splitCsv(header: true)
}

/*
 * Helper: read a samplesheet CSV and return the list of sample names
 */
def readSampleNames(samplesheet_path) {
    return readSamplesheet(samplesheet_path).collect { row -> row.sample_name }
}

/*
 * Helper: append an extension to a file path; returns a Path
 */
def appendExtension(path, extension) {
    return file(path).resolveSibling(file(path).name + extension)
}

/*
 * Helper: replace the file extension; returns a Path (e.g. test.fa -> test.dict)
 */
def replaceExtension(path, extension) {
    return file(path).resolveSibling(file(path).baseName + extension)
}

/*
 * Helper: resolve the pre-computed .paths ref-paths file from a Giraffe index
 * directory.  The file shares the same prefix as the .gbz.
 */
def resolveRefPathsFromIndexDir(idxDir) {
    def baseDir = file(idxDir, type: 'dir', checkIfExists: true)
    def gbzFiles = baseDir.list()?.findAll { it -> it.endsWith('.gbz') }?.sort()
    if (!gbzFiles) {
        error("No .gbz file found in Giraffe index directory: ${idxDir}")
    }
    def prefix = gbzFiles.first().replaceFirst(/\.gbz$/, '')
    return file(baseDir.resolve("${prefix}.ref_paths_dict"), checkIfExists: true)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow XOOS {
    take:
    ch_run_sheet // channel: run sheet read in from --input

    main:

    // =====================================================================
    // Parameter validation
    // =====================================================================

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Channel used to gather metric report JSON files from different steps across different
    // samples and runs to be used as input for the final metrics report generation step
    // during which we create one combined report for each run.
    // channel: [ val(run_meta), path(report) ]
    ch_run_metrics_report_json = channel.empty()

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
    // channel: [ val(meta), path(CollectGcBiasMetrics) ]
    ch_gc_bias_metrics_files = channel.empty()

    // Hybrid selection (HS) metrics files for publishing
    // channel: [ val(meta), path(*_metrics) ]
    ch_hs_metrics_files = channel.empty()

    // Read collapser metrics files for publishing (summary_stats.tsv,
    // cluster_size_distribution_summary.tsv, cluster_size_distributions.tsv)
    // channel: [ val(meta), path(read_collapser_metrics_dir) ]
    ch_read_collapser_metrics_files = channel.empty()

    // =====================================================================
    // STEP 1: Split input runs by file_type into appropriate channels
    // =====================================================================

    // Validate that run_dir and samplesheet exist before any processing.
    // Convert sectioned samplesheets ([Samples] format) to canonical flat CSV.
    ch_run_sheet = ch_run_sheet.map { run_meta, run_dir, samplesheet ->
        def validated_samplesheet = file(samplesheet, checkIfExists: true)
        def converted_samplesheet = convertSamplesheet(validated_samplesheet)
        validateSamplesheet(converted_samplesheet)
        [run_meta, file(run_dir, type: 'dir', checkIfExists: true), converted_samplesheet]
    }

    // Build a channel mapping each run to its list of sample names (parsed from the samplesheet).
    // This is used to set meta.sample_name on run-level (global) tasks so that the nf-report
    // plugin can attribute failures of these tasks to all affected samples.
    ch_run_sample_names = ch_run_sheet.map { run_meta, _run_dir, samplesheet ->
        [run_meta.run_name, readSampleNames(samplesheet)]
    }

    // Filter basecalls, inject 'id' and 'sample_name' list, and route to the subsampler.
    // Reuse ch_run_sample_names (computed once above) via combine to avoid re-parsing the samplesheet.
    ch_basecall_rdb = ch_run_sheet
        .filter { run_meta, _run_dir, _samplesheet -> run_meta.file_type == 'basecall_rdb' }
        .map { run_meta, run_dir, sheet -> [run_meta.run_name, run_meta + [id: run_meta.run_name], run_dir, sheet] }
        .combine(ch_run_sample_names, by: 0)
        .map { _run_name, run_meta, run_dir, sheet, sample_names -> [run_meta + [sample_name: sample_names], run_dir, sheet] }

    ch_basecall_fastq = ch_run_sheet
        .filter { run_meta, _run_dir, _samplesheet -> run_meta.file_type == 'basecall_fastq' }
        .map { run_meta, run_dir, sheet -> [run_meta.run_name, run_meta + [id: run_meta.run_name], run_dir, sheet] }
        .combine(ch_run_sample_names, by: 0)
        .map { _run_name, run_meta, run_dir, sheet, sample_names -> [run_meta + [sample_name: sample_names], run_dir, sheet] }

    // Run the subworkflow
    SUBSAMPLE_BASECALL_RUNS(
        ch_basecall_rdb,
        ch_basecall_fastq
    )

    // Remove the 'file_type' key before passing to demultiplexing, matching existing logic
    ch_basecall_files = SUBSAMPLE_BASECALL_RUNS.out.basecalls
        .map { run_meta, run_dir, samplesheet -> [removeKey(run_meta, 'file_type'), run_dir, samplesheet] }

    // Pre-demultiplexed FASTQ directories -> join with samplesheet metadata
    ch_demuxed_fastq_dirs = ch_run_sheet
        .filter { run_meta, _run_dir, _samplesheet -> run_meta.file_type == 'demultiplexed_fastq' }
        .map { run_meta, run_dir, samplesheet -> [removeKey(run_meta, 'file_type'), run_dir, samplesheet] }

    ch_demuxed_sample_meta = ch_demuxed_fastq_dirs.flatMap { run_meta, _run_dir, samplesheet ->
        readSamplesheet(samplesheet).collect { sample_meta -> [run_meta, sample_meta] }
    }
    ch_demultiplexed_fastqs = ch_demuxed_fastq_dirs
        .combine(ch_demuxed_sample_meta, by: 0)
        .map { run_meta, run_dir, _samplesheet, sample_meta ->
            def sample_fastq_files = []
            // Recursively search for FASTQ files in the demultiplexed FASTQ directory that match the sample name or sample SID in the metadata
            file(run_dir).eachFileRecurse { f ->
                if (f.name.endsWith('.fastq.gz') && (f.name.contains(sample_meta.sample_name) || f.name.contains(sample_meta.sample_sid))) {
                    sample_fastq_files << file(f, checkIfExists: true)
                }
            }
            [run_meta + sample_meta + [id: sample_meta.sample_name], sample_fastq_files]
        }
        .filter { _meta, sample_fastq_files -> sample_fastq_files.size() > 0 }

    // Pre-aligned BAM files -> skip demux and alignment
    ch_bam_runs = ch_run_sheet
        .filter { run_meta, _run_dir, _samplesheet -> run_meta.file_type == 'bam' }
        .map { run_meta, run_dir, samplesheet -> [removeKey(run_meta, 'file_type'), run_dir, samplesheet] }

    ch_bam_sample_meta = ch_bam_runs.flatMap { run_meta, _run_dir, samplesheet ->
        readSamplesheet(samplesheet).collect { sample_meta -> [run_meta, sample_meta] }
    }
    ch_bams = ch_bam_runs
        .combine(ch_bam_sample_meta, by: 0)
        .map { run_meta, run_dir, _samplesheet, sample_meta ->
            def bam_files_matching_sample_name_or_sid = []
            // Recursively search for BAM files in the run directory that match the sample name or sample SID in the metadata
            file(run_dir).eachFileRecurse { f ->
                if (f.name.endsWith('.bam') && (f.name.contains(sample_meta.sample_name) || f.name.contains(sample_meta.sample_sid))) {
                    bam_files_matching_sample_name_or_sid << f
                }
            }
            [run_meta + sample_meta + [id: sample_meta.sample_name], bam_files_matching_sample_name_or_sid]
        }
        .filter { _meta, bam_files -> bam_files.size() > 0 }
    // Multiple BAM files per sample cause file name collisions downstream
    // (BCFTOOLS_CONCAT, metrics report concat). Collect every match across all
    // samples first so we can report all duplicates in a single error rather
    // than failing on the first one — saves the operator a re-run loop when
    // several samples in the same input directory have duplicates.
    ch_bams = ch_bams
        .toList()
        .flatMap { entries ->
            def duplicates = entries.findAll { entry -> entry[1].size() > 1 }
            if (duplicates) {
                def details = duplicates.collect { entry ->
                    def meta = entry[0]
                    def bam_files = entry[1]
                    "  - ${meta.sample_name} (sample_sid: ${meta.sample_sid}, run: ${meta.run_name}): ${bam_files.size()} BAM files: ${bam_files.join(', ')}"
                }.join('\n')
                error("Found ${duplicates.size()} sample(s) with multiple matching BAM files. Only one BAM file per sample is allowed; remove or rename the extras before re-running:\n${details}")
            }
            return entries
        }
        .flatMap { meta, bam_files ->
            // Decode the single bam_type enum (declared on the run sheet or
            // per-sample samplesheet) into the two internal meta fields that the
            // dedup/rescue gates consume: bam_dedup and bam_secondary_rescued.
            // Only the rescue variants need unpacking; every other value maps
            // straight onto bam_dedup. Defaults to 'sas_markdup' for external
            // AXELIOS inputs that do not declare a state. This is the only place
            // the combined enum is interpreted — downstream steps read and write
            // the two fields directly via plain map merges.
            def bam_type = meta.bam_type ?: 'sas_markdup'
            def bam_secondary_rescued = bam_type.startsWith('xoos_secondary_rescue_and_')
            // Map each bam_type enum value to its underlying bam_dedup state.
            // This must stay in sync with the bam_type enum in
            // assets/schema_input.json and assets/schema_samplesheet.json. The
            // lookup fails fast on any value the schema did not catch rather than
            // silently passing an unrecognized string through as bam_dedup.
            def bam_dedup_by_type = [
                'no_dedup'                            : 'none',
                'sas_markdup'                         : 'sas_markdup',
                'pb_markdup'                          : 'pb_markdup',
                'xoos_markdup'                        : 'xoos_markdup',
                'xoos_consensus'                      : 'xoos_consensus',
                'xoos_secondary_rescue_and_markdup'   : 'xoos_markdup',
            ]
            def bam_dedup = bam_dedup_by_type[bam_type]
            if (bam_dedup == null) {
                error("Unknown bam_type '${bam_type}' for sample '${meta.sample_name}'. Allowed values are defined by the bam_type enum in assets/schema_input.json.")
            }
            def decoded_meta = removeKey(meta, 'bam_type') + [bam_dedup: bam_dedup, bam_secondary_rescued: bam_secondary_rescued]
            // Flatten list of BAM files so that each BAM file is a separate entry in the channel, and add corresponding BAI file -> [meta, bam, bai]
            bam_files.collect { bam_file ->
                [decoded_meta, file(bam_file, checkIfExists: true), file(appendExtension(bam_file, '.bai'), checkIfExists: true)]
            }
        }

    // Optionally reheader BAM-start inputs so their @RG SM/ID tags match the
    // samplesheet sample_name (meta.id). Only applies to pre-aligned BAMs;
    // freshly aligned BAMs already carry a correct @RG stamped during
    // alignment. SAMTOOLS_REHEADER emits the BAM only, so re-index with
    // SAMTOOLS_INDEX and rejoin to restore the [meta, bam, bai] shape that
    // downstream steps consume.
    if (params.enable_bam_reheader.toBoolean()) {
        SAMTOOLS_REHEADER(ch_bams)
        SAMTOOLS_INDEX_REHEADER(SAMTOOLS_REHEADER.out.bam)
        // join by meta yields [meta, bam, bai] directly
        ch_bams = SAMTOOLS_REHEADER.out.bam.join(SAMTOOLS_INDEX_REHEADER.out.index, by: 0)
    }

    // =====================================================================
    // STEP 2: Demultiplex basecall files
    // =====================================================================

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

    // =====================================================================
    // STEP 3: Alignment
    // =====================================================================

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

    // Build reference channels with [meta, path...] tuple structure
    ch_reference = params.reference_fasta
        ? channel.value([[id: 'reference'], file(params.reference_fasta, checkIfExists: true), file(params.reference_fai, checkIfExists: true)])
        : channel.value([[], [], []])

    ch_bwa_idx = params.aligner in ['pb_fq2bam', 'bwa']
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

    // --ref-paths for pbrun giraffe so the BAM @SQ header order matches the
    // reference used by downstream variant callers.
    // Without personalization: use the pre-built .dict derived from the linearized reference.
    // With personalization: use the pre-computed .paths file from the
    // unfiltered giraffe index directory (same prefix as the .gbz).
    if (params.aligner == 'pb_giraffe' && params.enable_pangenome_personalization.toBoolean()) {
        ch_giraffe_ref_paths = channel.value(resolveRefPathsFromIndexDir(params.unfiltered_giraffe_index_dir))
    } else if (params.aligner == 'pb_giraffe' && params.giraffe_reference_linearized) {
        ch_giraffe_ref_paths = channel.value(file(replaceExtension(params.giraffe_reference_linearized, '.ref_paths_dict'), checkIfExists: true))
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
    // Stamp the dedup state of freshly aligned BAMs. Parabricks aligners
    // (pb_fq2bam/pb_minimap2/pb_giraffe) mark duplicates during alignment, so
    // their output is 'pb_markdup' — except when alignment partitioning is
    // active, in which case per-partition markdup is disabled and the merged BAM
    // is unmarked ('none'), matching the markdup-on-merge behaviour in dedup.
    // All other aligners produce raw, unmarked BAMs ('none'). bam_secondary_rescued
    // is always false here since rescue only runs inside the dedup subworkflow.
    def aligner_dedup = (params.aligner in ['pb_fq2bam', 'pb_minimap2', 'pb_giraffe'] && !partitioningActive()) ? 'pb_markdup' : 'none'
    ch_bams = ch_bams.mix(
        align_out.bams.map { meta, bam, bai -> [removeKey(meta, 'personalized_on') + [bam_dedup: aligner_dedup, bam_secondary_rescued: false], bam, bai] }
    )

    // =====================================================================
    // STEP 4: Build shared resource channels for metrics and variant calling
    // =====================================================================

    ch_autosomes_bed = params.autosomes_bed
        ? channel.value([[id: 'autosomes_bed'], file(params.autosomes_bed, checkIfExists: true)])
        : channel.value([[], []])

    ch_target_bed = params.target_bed
        ? channel.value([[id: 'target_bed'], file(params.target_bed, checkIfExists: true)])
        : channel.value([[], []])

    ch_capture_bed = params.capture_bed
        ? channel.value([[id: 'capture_bed'], file(params.capture_bed, checkIfExists: true)])
        : channel.value([[id: 'capture_bed'], []])

    ch_high_confidence_bed_dir = params.high_confidence_bed_dir
        ? channel.value([[id: 'high_confidence_bed_dir'], file(params.high_confidence_bed_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], [], [], []])

    ch_fallback_high_confidence_bed = params.fallback_high_confidence_bed
        ? channel.value([[id: 'fallback_high_confidence_bed'], file(params.fallback_high_confidence_bed, checkIfExists: true)])
        : channel.value([[], []])

    ch_tfe_included_regions = params.tfe_included_regions_bed
        ? channel.value([[id: 'tfe_included_regions'], file(params.tfe_included_regions_bed, checkIfExists: true)])
        : channel.value([[], []])

    ch_pop_af_vcf_with_tbi = params.pop_af_vcf
        ? channel.value([[id: 'pop_af_vcf'], file(params.pop_af_vcf, checkIfExists: true), file(params.pop_af_tbi, checkIfExists: true)])
        : channel.value([[], [], []])

    ch_contamination_vcf = params.contamination_vcf
        ? channel.value([[id: 'contamination_vcf'], file(params.contamination_vcf, checkIfExists: true)])
        : channel.value([[], []])

    ch_contamination_included_regions = params.contamination_included_regions_bed
        ? channel.value([[id: 'contamination_included_regions'], file(params.contamination_included_regions_bed, checkIfExists: true)])
        : channel.value([[], []])

    // Picard reference: use linearized ref for giraffe aligners, otherwise the main reference
    def candidate_linearized_ref = params.enable_pangenome_personalization.toBoolean()
        ? params.unfiltered_giraffe_reference_linearized
        : params.giraffe_reference_linearized
    def linearized_ref = params.aligner in ['pb_giraffe', 'vg_giraffe'] && candidate_linearized_ref
        ? candidate_linearized_ref
        : params.reference_fasta
    ch_picard_reference = linearized_ref
        ? channel.value([[id: 'picard_reference'], file(linearized_ref, checkIfExists: true), file(appendExtension(linearized_ref, '.fai'), checkIfExists: true)])
        : channel.value([[], [], []])
    // Policy for secondary-alignment rescue: whether this run's configuration
    // wants rescue at all. Whether it actually runs for a given BAM is decided
    // per-sample inside DEDUP, which skips BAMs already carrying
    // bam_secondary_rescued (rescue is not idempotent — re-running it on an
    // already-rescued BAM corrupts the alignment record set).
    def rescue_policy = (params.aligner in ['vg_giraffe', 'pb_giraffe'] && params.enable_copy_number_caller.toBoolean())

    // =====================================================================
    // STEP 5: Pre-dedup metrics (optional)
    // =====================================================================

    if (params.enable_pre_dedup_metrics.toBoolean() || (params.dedup_strategy == 'none' && (!params.enable_subsample_bam.toBoolean() || params.enable_pre_subsample_metrics.toBoolean()))) {
        // Calculate pre-dedup metrics if
        // - dedup is diabled and either subsampling is disabled or pre-subsample metrics are explicitly enabled, OR
        // - dedup is enabled and pre-dedup metrics are explicitly enabled
        // If both dedup and subsampling are disabled, the pre-dedup BAM is the final BAM used for downstream analysis
        if (params.aligner in ['pb_fq2bam', 'pb_minimap2', 'pb_giraffe'] && params.dedup_strategy == 'markdup' && !partitioningActive() && !params.enable_read_collapser_dedup.toBoolean()) {
            log.warn("Pre-dedup metrics cannot be generated when using Parabricks aligners since duplicate marking is integrated into the alignment process.")
        }
        else {
            PRE_DEDUP_METRICS(
                'original',
                ch_bams,
                ch_reference,
                ch_autosomes_bed,
                ch_target_bed,
                ch_capture_bed,
                ch_high_confidence_bed_dir,
                ch_fallback_high_confidence_bed,
                ch_picard_reference,
                params.enable_alignment_metrics_autosomes.toBoolean(),
                params.enable_alignment_metrics_targets.toBoolean(),
                params.enable_alignment_metrics_high_confidence.toBoolean(),
                params.enable_gc_bias_metrics.toBoolean(),
                params.enable_hs_metrics.toBoolean(),
            )
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(PRE_DEDUP_METRICS.out.report)
        }
    }

    // =====================================================================
    // STEP 6: Dedup (markdup or consensus)
    // =====================================================================

    if (params.dedup_strategy != 'none') {
        DEDUP(
            params.dedup_strategy,
            params.aligner,
            params.enable_read_collapser_dedup.toBoolean(),
            rescue_policy,
            ch_bams,
            ch_target_bed,
            ch_reference,
            ch_bwa_idx,
            ch_giraffe_idx,
            ch_giraffe_ref_paths,
        )
        ch_bams = DEDUP.out.bam
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
        if (params.enable_pre_subsample_metrics.toBoolean() || !params.enable_subsample_bam.toBoolean()) {
            // Calculate pre-subsample dedup metrics if subsampling is disabled or if pre-subsample metrics are explicitly enabled
            // If subsampling is disabled, the post-dedup BAM is the final BAM used for downstream analysis
            DEDUP_METRICS(
                'dedup',
                ch_bams,
                ch_reference,
                ch_autosomes_bed,
                ch_target_bed,
                ch_capture_bed,
                ch_high_confidence_bed_dir,
                ch_fallback_high_confidence_bed,
                ch_picard_reference,
                params.enable_alignment_metrics_autosomes.toBoolean(),
                params.enable_alignment_metrics_targets.toBoolean(),
                params.enable_alignment_metrics_high_confidence.toBoolean(),
                params.enable_gc_bias_metrics.toBoolean(),
                params.enable_hs_metrics.toBoolean(),
            )
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(DEDUP_METRICS.out.report)
        }
    }

    // =====================================================================
    // STEP 7: Subsample BAM by coverage (optional)
    // =====================================================================

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
        // as this is the final output BAM that will be used for downstream analysis
        SUBSAMPLE_METRICS(
            params.dedup_strategy != 'none' ? 'dedup_subsample' : 'subsample',
            ch_bams,
            ch_reference,
            ch_autosomes_bed,
            ch_target_bed,
            ch_capture_bed,
            ch_high_confidence_bed_dir,
            ch_fallback_high_confidence_bed,
            ch_picard_reference,
            params.enable_alignment_metrics_autosomes.toBoolean(),
            params.enable_alignment_metrics_targets.toBoolean(),
            params.enable_alignment_metrics_high_confidence.toBoolean(),
            params.enable_gc_bias_metrics.toBoolean(),
            params.enable_hs_metrics.toBoolean(),
        )
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(SUBSAMPLE_METRICS.out.report)
    }

    // =====================================================================
    // STEP 7.5: Contamination Detection
    // =====================================================================

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

    // =====================================================================
    // STEP 7.75: Tumor Fraction Estimation
    // =====================================================================

    ch_tfe_output = channel.empty()

    if (params.enable_tfe.toBoolean()) {
        TUMOR_FRACTION_ESTIMATION(
            ch_bams,
            ch_reference,
            ch_tfe_included_regions,
        )
        ch_tfe_output = TUMOR_FRACTION_ESTIMATION.out.results
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
            XOOS_PARSE_TUMOR_FRACTION_METRICS(TUMOR_FRACTION_ESTIMATION.out.results).report
        )
    }

    // =====================================================================
    // STEP 8: Variant calling
    // =====================================================================

    ch_vcf = channel.empty()
    ch_raw_vc_vcf = channel.empty()
    ch_svc_pretrained_vcf = channel.empty()
    ch_svc_retrained_vcf = channel.empty()
    ch_vcfeval_output = channel.empty()
    ch_cnc_output = channel.empty()
    ch_str_output = channel.empty()

    // Build reference resource channels for variant calling.
    // Giraffe aligners produce BAMs against the pangenome linearized reference
    // (179 contigs) rather than the standard no-alt reference (2580 contigs).
    // Use the matching linearized giraffe reference for variant calling so the
    // reference dictionary is consistent with the BAM. This mirrors the pattern
    // used for Picard metrics (ch_picard_reference).
    def candidate_linearized_vc_ref = params.enable_pangenome_personalization.toBoolean()
        ? params.unfiltered_giraffe_reference_linearized
        : params.giraffe_reference_linearized
    def use_linearized_vc_ref = params.aligner in ['pb_giraffe', 'vg_giraffe'] && candidate_linearized_vc_ref
    ch_vc_reference = use_linearized_vc_ref
        ? channel.value([[id: 'vc_reference'], file(candidate_linearized_vc_ref, checkIfExists: true), file(appendExtension(candidate_linearized_vc_ref, '.fai'), checkIfExists: true)])
        : ch_reference
    ch_vc_reference_dict = use_linearized_vc_ref
        ? channel.value([[id: 'vc_reference_dict'], file(replaceExtension(candidate_linearized_vc_ref, '.dict'), checkIfExists: true)])
        : (params.reference_dict
            ? channel.value([[id: 'vc_reference_dict'], file(params.reference_dict, checkIfExists: true)])
            : channel.value([[], []]))

    // Subset reference used ONLY for variant-calling partitioning (GATK4_SPLITINTERVALS
    // and the FAI_TO_BED whole-genome interval fallback). It is graph-specific: when
    // aligning to a graph (personalized or not) we already know the base-graph folder
    // via candidate_linearized_vc_ref. By convention the linearized reference is
    // <prefix>.fa and the partitioning subset is <prefix>.subset.fa in the same folder,
    // with sidecars <prefix>.subset.fa.fai and <prefix>.subset.dict. The subset contains
    // a smaller set of chromosomes so partitioning skips _random and other contigs. The
    // variant callers themselves still use the full ch_vc_reference. When the aligner is
    // non-graph or no subset file is present, fall back to the full reference so existing
    // behavior is unchanged.
    // Derive the subset path from the reference URI string (replacing the trailing
    // .fa with .subset.fa) and resolve it with file(checkIfExists: true). This uses the
    // same mechanism Nextflow uses for the .fai/.dict sidecars above, which is reliable
    // for remote references (e.g. s3://); Groovy's .exists() on a resolveSibling Path can
    // silently report false for remote paths. A missing subset is not an error here, so
    // the lookup is wrapped to fall back to the full reference.
    def candidate_subset_str = use_linearized_vc_ref
        ? candidate_linearized_vc_ref.toString().replaceAll(/\.fa$/, '.subset.fa')
        : null
    def candidate_subset_ref = null
    if (use_linearized_vc_ref) {
        try {
            candidate_subset_ref = file(candidate_subset_str, checkIfExists: true)
        }
        catch (Exception e) {
            candidate_subset_ref = null
        }
    }
    def use_subset_ref = candidate_subset_ref != null
    if (use_linearized_vc_ref) {
        if (use_subset_ref) {
            log.info("Variant-calling partitioning will use subset reference: ${candidate_subset_ref}")
        }
        else {
            log.warn("Graph aligner '${params.aligner}' is in use but no subset reference was found at ${candidate_subset_str}; partitioning will fall back to the full reference. To restrict partitioning to a smaller chromosome set, place <prefix>.subset.fa (with .subset.fa.fai and .subset.dict sidecars) next to the linearized reference.")
        }
    }
    ch_splitintervals_reference = use_subset_ref
        ? channel.value([[id: 'splitintervals_reference'], candidate_subset_ref, file(appendExtension(candidate_subset_ref, '.fai'), checkIfExists: true)])
        : ch_vc_reference
    ch_splitintervals_reference_dict = use_subset_ref
        ? channel.value([[id: 'splitintervals_reference_dict'], file(replaceExtension(candidate_subset_ref, '.dict'), checkIfExists: true)])
        : ch_vc_reference_dict

    ch_germline_svc_truth_vcf_dir = params.germline_small_variant_truth_vcf_dir
        ? channel.value([[id: 'truth_vcf_dir'], file(params.germline_small_variant_truth_vcf_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])

    ch_germline_svc_train_bed_dir = params.germline_small_variant_train_bed_dir
        ? channel.value([[id: 'train_bed_dir'], file(params.germline_small_variant_train_bed_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])

    ch_germline_svc_test_bed_dir = params.germline_small_variant_test_bed_dir
        ? channel.value([[id: 'test_bed_dir'], file(params.germline_small_variant_test_bed_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])

    ch_germline_svc_eval_bed_dir = params.germline_small_variant_eval_bed_dir
        ? channel.value([[id: 'eval_bed_dir'], file(params.germline_small_variant_eval_bed_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])

    ch_somatic_tn_svc_truth_vcf_dir = params.somatic_tumor_normal_small_variant_truth_vcf_dir
        ? channel.value([[id: 'somatic_tn_truth_vcf_dir'], file(params.somatic_tumor_normal_small_variant_truth_vcf_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])

    ch_somatic_tn_svc_eval_bed_dir = params.somatic_tumor_normal_small_variant_eval_bed_dir
        ? channel.value([[id: 'somatic_tn_eval_bed_dir'], file(params.somatic_tumor_normal_small_variant_eval_bed_dir, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])

    ch_svc_germline_models = params.svc_model_germline_snv && params.svc_model_germline_indel
        ? channel.value([[id: 'svc_models_germline'], file(params.svc_model_germline_snv, checkIfExists: true), file(params.svc_model_germline_indel, checkIfExists: true)])
        : channel.value([[], [], []])

    ch_svc_somatic_tn_models = params.svc_model_somatic_tumor_normal
        ? channel.value([[id: 'svc_model_somatic_tumor_normal'], file(params.svc_model_somatic_tumor_normal, checkIfExists: true)])
        : channel.value([[], []])

    ch_par_beds = params.chr_x_par_bed
        ? channel.value([[id: 'chr_x_par_bed'], file(params.chr_x_par_bed, checkIfExists: true), file(params.chr_y_par_bed, checkIfExists: true)])
        : channel.value([[], []])

    ch_mappability_bw = params.mappability_bw
        ? channel.value([[id: 'mappability_bw'], file(params.mappability_bw, checkIfExists: true)])
        : channel.value([[], []])

    ch_cnv_seed_file = params.cnv_seed_seg
        ? channel.value([[id: 'cnv_seed_file'], file(params.cnv_seed_seg, checkIfExists: true)])
        : channel.value([[], []])

    ch_cnv_blocklist_bed = params.cnv_blocklist_bed
        ? channel.value([[id: 'cnv_blocklist_bed'], file(params.cnv_blocklist_bed, checkIfExists: true)])
        : channel.value([[], []])

    ch_repeat_catalog = params.repeat_catalog
        ? channel.value([[id: 'repeat_catalog'], file(params.repeat_catalog, checkIfExists: true)])
        : channel.value([[], []])

    def reference_sdf = params.reference_sdf ?: null

    // Enable variant calling partitioning if the option is enabled and the partition count for the
    // chosen variant caller is greater than 1.
    def use_cpu_variant_callers =
        (params.variant_calling_mode == 'somatic_tn' && params.use_gatk_mutect2.toBoolean()) ||
        (params.variant_calling_mode == 'germline' && params.use_gatk_haplotype_caller.toBoolean())
    def vc_partition_count = use_cpu_variant_callers
        ? params.cpu_variant_calling_partition_count.toInteger()
        : params.gpu_variant_calling_partition_count.toInteger()
    def vc_partitioning_enabled = vc_partition_count > 1 && params.enable_variant_calling_partitioning.toBoolean()

    if (params.variant_calling_mode == 'somatic_tn') {
        ch_samples_by_group = ch_bams
            .filter { meta, _bam, _bai -> meta.containsKey('TN') && meta.containsKey('group') && meta.TN && meta.group }
            // error() is intentional here — do not downgrade to log.warn().
            // ifEmpty injects the closure's return value into the channel when empty;
            // error() throws, preventing a null from reaching the downstream .map/.groupTuple.
            .ifEmpty { error("variant_calling_mode is 'somatic_tn' but no samples have 'TN' and 'group' columns set in the samplesheet. Somatic tumor-normal variant calling requires paired samples with TN='T'/'N' and a shared group identifier.") }
            // The pairing key is deliberately reduced to [group:] only — it must
            // NOT use the full meta. Tumor and normal BAMs can legitimately carry
            // different per-sample dedup state (e.g. bam_dedup/bam_secondary_rescued),
            // so grouping on full meta would split a pair into two singletons.
            .map { meta, bam, bai -> [[group: meta.group], [meta, bam, bai]] }
            .groupTuple(size: 2, remainder: false)
        // Group BAM files into pairs and discard any unpaired BAM files
        ch_tn_bam_pairs = ch_samples_by_group
            .filter { _group, samples -> samples.size() == 2 }
            .map { group, samples ->
                def tumor = samples.find { meta, _bam, _bai -> meta.TN && meta.TN.toString().toUpperCase() == 'T' }
                def normal = samples.find { meta, _bam, _bai -> meta.TN && meta.TN.toString().toUpperCase() == 'N' }
                if (tumor && normal) {
                    def meta = tumor[0].clone()
                    meta.remove('TN')
                    meta.remove('sample_name')
                    meta.remove('sample_sid')
                    meta.remove('group')
                    // Drop per-sample dedup state when collapsing the pair into a
                    // single meta. These keys describe an individual BAM's
                    // dedup/rescue history, which is meaningless at pair level and
                    // could otherwise propagate only the tumor's state to
                    // downstream somatic logic. Dedup is already complete here.
                    meta.remove('bam_dedup')
                    meta.remove('bam_secondary_rescued')
                    meta.id = group.group
                    meta.tumor_id = tumor[0].id
                    meta.normal_id = normal[0].id

                    // Resolve sample_type for model selection (default: ffpe).
                    // If both tumor and normal specify sample_type but disagree, warn and fall back to ffpe.
                    def tumor_sample_type = tumor[0].sample_type ?: 'ffpe'
                    def normal_sample_type = normal[0].sample_type ?: 'ffpe'
                    if (tumor_sample_type != normal_sample_type) {
                        log.warn("sample_type mismatch for group ${group.group}: tumor='${tumor_sample_type}', normal='${normal_sample_type}'. Defaulting to 'ffpe'.")
                        meta.sample_type = 'ffpe'
                    } else {
                        meta.sample_type = tumor_sample_type
                    }

                    // Return [meta, [tumor_bam, normal_bam], [tumor_bai, normal_bai]]
                    return [meta, [tumor[1], normal[1]], [tumor[2], normal[2]]]
                }
                else {
                    log.warn("Could not identify tumor and normal BAM files for group ${group.group}.")
                    null
                }
            }
            .filter { it -> it != null }

        CALL_SOMATIC_VARIANTS_TN(
            ch_tn_bam_pairs,
            ch_vc_reference,
            ch_vc_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_target_bed,
            ch_pop_af_vcf_with_tbi,
            ch_svc_somatic_tn_models,
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
            params.aligner,
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
        ch_svc_pretrained_vcf = CALL_SOMATIC_VARIANTS_TN.out.svc_vcf
        ch_vcfeval_output = CALL_SOMATIC_VARIANTS_TN.out.small_variant_vcfeval_output
        ch_cnc_output = CALL_SOMATIC_VARIANTS_TN.out.copy_number_caller_output
        ch_run_metrics_report_json = ch_run_metrics_report_json.mix(CALL_SOMATIC_VARIANTS_TN.out.small_variant_vcfeval_report)

        if (params.enable_copy_number_caller.toBoolean()) {
            ch_run_metrics_report_json = ch_run_metrics_report_json.mix(
                XOOS_PARSE_CNC_METRICS_SOMATIC(CALL_SOMATIC_VARIANTS_TN.out.copy_number_caller_output).report
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
            ch_repeat_catalog,
            params.enable_small_variant_caller.toBoolean(),
            params.use_gatk_haplotype_caller.toBoolean(),
            params.use_multisample_model.toBoolean(),
            params.aligner,
            params.enable_small_variant_caller_retrain.toBoolean(),
            params.enable_vcfeval.toBoolean(),
            params.enable_copy_number_caller.toBoolean(),
            params.enable_copy_number_caller_baf.toBoolean(),
            params.enable_str_caller.toBoolean(),
            params.enable_merge_vcf.toBoolean(),
            reference_sdf,
            vc_partitioning_enabled,
            vc_partition_count,
        )
        ch_vcf = CALL_GERMLINE_VARIANTS.out.merged_vcf
        ch_raw_vc_vcf = CALL_GERMLINE_VARIANTS.out.raw_vc_vcf
            .map { meta, vcf, tbi -> [meta + [raw_variant_caller: 'haplotype_caller'], vcf, tbi] }
        ch_svc_pretrained_vcf = CALL_GERMLINE_VARIANTS.out.svc_pretrained_vcf
        ch_svc_retrained_vcf = CALL_GERMLINE_VARIANTS.out.svc_retrained_vcf
        ch_vcfeval_output = CALL_GERMLINE_VARIANTS.out.small_variant_vcfeval_output
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

    // =====================================================================
    // STEP 9: Metrics reporting
    // =====================================================================

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
    else if (params.dedup_strategy != 'none') {
        ch_custom_multiqc_files_with_meta = ch_custom_multiqc_files_with_meta.mix(DEDUP_METRICS.out.custom_multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_multiqc_files_with_meta = ch_multiqc_files_with_meta.mix(DEDUP_METRICS.out.multiqc_files.map { meta, metrics_dir -> [[id: meta.run_name, run_type: meta.run_type], metrics_dir] })
        ch_alignment_metrics_files = ch_alignment_metrics_files.mix(DEDUP_METRICS.out.alignment_metrics_files)
        ch_gc_bias_metrics_files = ch_gc_bias_metrics_files.mix(DEDUP_METRICS.out.gc_bias_metrics_files)
        ch_hs_metrics_files = ch_hs_metrics_files.mix(DEDUP_METRICS.out.hs_metrics_files)
    }
    else {
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
    pre_subsample_bam       = ch_pre_subsample_bam // channel: [ [meta], path(bam), path(bai) ]
    vcf                     = ch_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    raw_vc_vcf              = ch_raw_vc_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    svc_pretrained_vcf      = ch_svc_pretrained_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    svc_retrained_vcf       = ch_svc_retrained_vcf // channel: [ [meta], path(vcf), path(tbi) ]
    vcfeval_output          = ch_vcfeval_output // channel: [ [meta], path(vcfeval_dir) ]
    demuxed_fastqs          = ch_demuxed_fastqs_output // channel: [ [meta], path(sample_fastq_dir) ]
    demux_metrics_files     = ch_demux_metrics_files // channel: [ [meta], path(metrics_dir) ]
    cnc_output              = ch_cnc_output // channel: [ [meta], path(copy_number_caller_output) ]
    str_output              = ch_str_output // channel: [ [meta], path(str_caller_output) ]
    tfe_output              = ch_tfe_output // channel: [ [meta], path(tfe_results_dir) ]
    contamination_output    = ch_contamination_output // channel: [ [meta], path(contamination_results_dir) ]
    alignment_metrics_files = ch_alignment_metrics_files // channel: [ [meta], path(read_metrics), path(coverage), path(accuracy) ]
    gc_bias_metrics_files   = ch_gc_bias_metrics_files // channel: [ [meta], path(CollectGcBiasMetrics) ]
    hs_metrics_files        = ch_hs_metrics_files // channel: [ [meta], path(*_metrics) ]
    read_collapser_metrics_files = ch_read_collapser_metrics_files // channel: [ [meta], path(read_collapser_metrics_dir) ]
    combined_report         = ch_combined_metrics_report // channel: [ val(run_meta), path(report) ]
    combined_report_summary = ch_combined_metrics_report_summary // channel: [ val(run_meta), path(report_summary) ]
}
