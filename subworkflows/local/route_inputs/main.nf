include { SUBSAMPLE_BASECALL_RUNS        } from '../subsample_basecall_runs/main'
include { SAMTOOLS_REHEADER              } from '../../../modules/nf-core/samtools/reheader/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_REHEADER } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CONVERT_CRAM_TO_BAM } from '../../../modules/nf-core/samtools/convert/main'
include { VERIFY_CRAM                    } from '../../../modules/local/verify_cram/main'
include { convertSamplesheet ; validateSamplesheet } from '../utils_samplesheet/main'
include { removeKey ; readSamplesheet ; readSampleNames ; appendExtension } from '../utils_xoos_workflow/main'

/*
 * Decode the bam_type enum into the internal metadata used by the dedup and
 * secondary-rescue gates. An omitted value deliberately defaults to no_dedup
 * so legacy inputs receive the configured deduplication strategy.
 */
def decodeBamType(meta) {
    def bam_type = meta.bam_type ?: 'no_dedup'
    def bam_dedup_type_by_bam_type = [
        'no_dedup'                          : 'none',
        'sas_markdup'                       : 'sas_markdup',
        'pb_markdup'                        : 'pb_markdup',
        'xoos_markdup'                      : 'xoos_markdup',
        'xoos_secondary_rescue_and_markdup' : 'xoos_markdup',
        'xoos_consensus'                    : 'xoos_consensus',
    ]
    def bam_dedup_type = bam_dedup_type_by_bam_type[bam_type]
    if (bam_dedup_type == null) {
        error("Unknown bam_type '${bam_type}' for sample '${meta.sample_name}'. Allowed values are defined by the bam_type enum in assets/schema_input.json.")
    }

    def required_dedup_strategy_by_bam_type = [
        'xoos_markdup'                      : 'markdup',
        'xoos_secondary_rescue_and_markdup' : 'markdup',
        'xoos_consensus'                    : 'consensus',
    ]
    def required_dedup_strategy = required_dedup_strategy_by_bam_type[bam_type]
    if (required_dedup_strategy && params.dedup_strategy != required_dedup_strategy) {
        if (params.dedup_strategy == 'consensus') {
            // Allow consensus to proceed with a warning even if the dedup_strategy does not match bam_type
            log.warn("Sample '${meta.sample_name}' has bam_type '${bam_type}' but dedup_strategy='${params.dedup_strategy}'. This often indicates a possible misconfiguration of dedup_strategy.")
        } else {
            // In all other cases (bam_type='xoos_markdup/xoos_secondary_rescue_and_markdup' and dedup_strategy='consensus/none')
            // fail because we cannot run markdup on a consensus BAM or run no dedup on a markdup BAM
            error("Sample '${meta.sample_name}' has bam_type '${bam_type}' which requires pipeline dedup_strategy='${required_dedup_strategy}', but dedup_strategy='${params.dedup_strategy}'.")
        }
    }

    return removeKey(meta, 'bam_type') + [
        bam_dedup_type: bam_dedup_type,
        bam_secondary_rescued: bam_type == 'xoos_secondary_rescue_and_markdup',
    ]
}

/*
 * Discover pre-aligned input files (BAM or CRAM) for each sample.
 *
 * Recursively searches each run directory for files with the given extension
 * that match the sample_name or sample_sid, fails fast when a sample matches
 * more than one file (which would collide downstream), decodes the bam_type
 * enum, and resolves the sidecar index. Shared by the BAM and CRAM routing
 * paths so both discover files with identical semantics.
 *
 * Returns a channel of [meta, file, index] tuples.
 */

def discoverAlignedInputs(ch_runs, String extension, String indexExtension) {
    def ch_runs_no_type = ch_runs
        .map { run_meta, run_dir, samplesheet -> [removeKey(run_meta, 'file_type'), run_dir, samplesheet] }

    def ch_sample_meta = ch_runs_no_type.flatMap { run_meta, _run_dir, samplesheet ->
        readSamplesheet(samplesheet).collect { sample_meta -> [run_meta, sample_meta] }
    }

    def ch_discovered = ch_runs_no_type
        .combine(ch_sample_meta, by: 0)
        .map { run_meta, run_dir, _samplesheet, sample_meta ->
            def files_matching_sample_name_or_sid = []
            // A non-empty per-sample value overrides the run-level state. An
            // empty optional samplesheet column must not erase the run default.
            def effective_sample_meta = sample_meta.bam_type ? sample_meta : removeKey(sample_meta, 'bam_type')
            // Recursively search the run directory for files that match the sample name or sample SID in the metadata
            file(run_dir).eachFileRecurse { f ->
                if (f.name.endsWith(extension) && (f.name.contains(sample_meta.sample_name) || f.name.contains(sample_meta.sample_sid))) {
                    files_matching_sample_name_or_sid << f
                }
            }
            [run_meta + effective_sample_meta + [id: sample_meta.sample_name], files_matching_sample_name_or_sid]
        }
        .filter { _meta, matched_files -> matched_files.size() > 0 }

    // Multiple files per sample cause file name collisions downstream
    // (BCFTOOLS_CONCAT, metrics report concat). Collect every match across all
    // samples first so we can report all duplicates in a single error rather
    // than failing on the first one — saves the operator a re-run loop when
    // several samples in the same input directory have duplicates.
    return ch_discovered
        .toList()
        .flatMap { entries ->
            def duplicates = entries.findAll { entry -> entry[1].size() > 1 }
            if (duplicates) {
                def details = duplicates.collect { entry ->
                    def meta = entry[0]
                    def matched_files = entry[1]
                    "  - ${meta.sample_name} (sample_sid: ${meta.sample_sid}, run: ${meta.run_name}): ${matched_files.size()} ${extension} files: ${matched_files.join(', ')}"
                }.join('\n')
                error("Found ${duplicates.size()} sample(s) with multiple matching ${extension} files. Only one ${extension} file per sample is allowed; remove or rename the extras before re-running:\n${details}")
            }
            return entries
        }
        .flatMap { meta, matched_files ->
            def decoded_meta = decodeBamType(meta)
            // Flatten list of files so that each file is a separate entry in the channel, and add corresponding index -> [meta, file, index]
            matched_files.collect { aligned_file ->
                [decoded_meta, file(aligned_file, checkIfExists: true), file(appendExtension(aligned_file, indexExtension), checkIfExists: true)]
            }
        }
}

/*
 * Route input runs by file_type into per-type channels.
 *
 * Handles five input types declared in the run sheet:
 *   basecall_rdb / basecall_fastq → subsampled then demultiplexed by DEMULTIPLEX
 *   demultiplexed_fastq           → FASTQ files that skip demux
 *   bam                           → pre-aligned BAMs that skip demux + alignment
 *   cram                          → pre-aligned CRAMs, verified to carry M5 tags
 *
 * Emits:
 *   basecall_files       – [meta, path(run_dir), path(samplesheet)]
 *   demultiplexed_fastqs – [meta, List<path(fastq)>]
 *   bams                 – [meta, path(bam), path(bai)]
 *   run_sample_names     – [run_name, List<sample_name>]
 */
workflow ROUTE_INPUTS {
    take:
    ch_run_sheet      // channel: [ val(run_meta), path(run_dir), path(samplesheet) ]
    ch_reference      // channel: [ val(meta), path(fasta), path(fai) ] (CRAM->BAM conversion)

    main:

    // Validate every per-run samplesheet on the main thread before any channel
    // operators run. samplesheetToList (inside readSamplesheet) throws
    // SchemaValidationException, which nf-schema only formats as a clean boxed
    // message when the call originates on the main thread. Calling it inside a
    // .map/.flatMap closure (actor thread) causes DataflowHelper to wrap the
    // exception in a raw reflection stack trace.
    file(params.input, checkIfExists: true).splitCsv(header: true).each { row ->
        if (row.samplesheet) {
            readSamplesheet(row.samplesheet.trim())
        }
    }

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

    // Pre-aligned BAM files -> skip demux and alignment; discover [meta, bam, bai].
    ch_bam_runs = ch_run_sheet
        .filter { run_meta, _run_dir, _samplesheet -> run_meta.file_type == 'bam' }
    ch_bams = discoverAlignedInputs(ch_bam_runs, '.bam', '.bai')

    // Pre-aligned CRAM files -> discover [meta, cram, crai], verify they carry
    // M5 tags, convert to BAM, and mix into ch_bams before the reheader step so
    // converted CRAMs are handled identically to pre-aligned BAM inputs.
    ch_cram_runs = ch_run_sheet
        .filter { run_meta, _run_dir, _samplesheet -> run_meta.file_type == 'cram' }
    ch_crams = discoverAlignedInputs(ch_cram_runs, '.cram', '.crai')
        .map { meta, cram, crai ->
            if (!params.reference_fasta || !params.reference_fai) {
                error("reference_fasta and reference_fai are required when file_type is 'cram'. CRAM files require the reference genome for conversion to BAM. Provide --reference_fasta and --reference_fai.")
            }
            [meta, cram, crai]
        }

    // CRAM-only steps are gated behind this flag to keep VERIFY_CRAM and
    // SAMTOOLS_CONVERT_CRAM_TO_BAM out of the log unless CRAM inputs are actually present.
    def has_cram_inputs = file(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .any { row -> row.file_type?.trim() == 'cram' }

    // Pre-aligned CRAM files -> discover [meta, cram, crai], verify they carry
    // M5 tags, convert to BAM, and mix into ch_bams before the reheader step so
    // converted CRAMs are handled identically to pre-aligned BAM inputs.
    if (has_cram_inputs) {
        ch_cram_runs = ch_run_sheet
            .filter { run_meta, _run_dir, _samplesheet -> run_meta.file_type == 'cram' }
        ch_crams = discoverAlignedInputs(ch_cram_runs, '.cram', '.crai')
            .map { meta, cram, crai ->
                if (!params.reference_fasta || !params.reference_fai) {
                    error("reference_fasta and reference_fai are required when file_type is 'cram'. CRAM files require the reference genome for conversion to BAM. Provide --reference_fasta and --reference_fai.")
                }
                [meta, cram, crai]
            }

        // Verify each pre-aligned CRAM carries M5 (reference-checksum) tags on
        // its @SQ header records so the reference used for encoding can be
        // resolved during conversion; fail fast rather than silently producing a
        // corrupt BAM.
        VERIFY_CRAM(ch_crams)

        SAMTOOLS_CONVERT_CRAM_TO_BAM(
            VERIFY_CRAM.out.cram,
            ch_reference,
        )
        ch_bams = ch_bams.mix(
            SAMTOOLS_CONVERT_CRAM_TO_BAM.out.bam.join(SAMTOOLS_CONVERT_CRAM_TO_BAM.out.bai)
        )
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

    emit:
    basecall_files        = ch_basecall_files        // channel: [ val(meta), path(run_dir), path(samplesheet) ]
    demultiplexed_fastqs  = ch_demultiplexed_fastqs  // channel: [ val(meta), List<path(fastq)> ]
    bams                  = ch_bams                  // channel: [ val(meta), path(bam), path(bai) ]
    run_sample_names      = ch_run_sample_names       // channel: [ val(run_name), List<val(sample_name)> ]
}
