include { XOOS_PARSE_DEMUX_METRICS      } from '../../../modules/local/xoos_generate_metrics_report/parse_demux_metrics/main'
include { XOOS_DEMUX_METRICS_TO_MULTIQC } from '../../../modules/local/xoos_metrics_to_multiqc/demux_metrics/main'
include { XOOS_DEMUX                    } from '../../../modules/local/xoos_demux/main'

def readSamplesheet(meta, samplesheet_path) {
    return file(samplesheet_path)
        .splitCsv(header: true)
        .collect { row -> [meta + [id: row.sample_name], row] }
}

workflow DEMULTIPLEX {
    take:
    ch_input // channel: [ val(run_level_meta), path(fastq_dir), path(samplesheet) ]
    ch_strand_detect_index // channel: path(strand_detect_index)

    main:
    demux_output = XOOS_DEMUX(
        ch_input,
        ch_strand_detect_index,
    )

    // Create a channel of demultiplexed samples along with their metadata
    // (run-level metadata from input channel + sample name from directory name)
    ch_samples = demux_output.demux_output
        .map { run_level_meta, dirs -> [run_level_meta, dirs.findAll { dir -> dir.isDirectory() && dir.name != 'metrics' }] }
        .flatMap { run_level_meta, sample_dirs -> sample_dirs.collect { sample_dir -> [run_level_meta + [id: sample_dir.name], sample_dir] } }
    // Combine run-level metadata with sample-level metadata from the samplesheet
    ch_sample_meta = ch_input.flatMap { run_level_meta, _fastq_dir, samplesheet_path -> readSamplesheet(run_level_meta, samplesheet_path) }
    // channel: [ val(run_level_meta_with_id), sample_meta ]
    ch_samples_with_meta = ch_samples
        .join(ch_sample_meta)
        .map { meta, sample_dir, additional_meta -> [meta + additional_meta, sample_dir] }
    // channel: [ val(combined_meta), path(sample_dir) ]

    // When exclude_partial_reads is enabled for simplex adapter types, demux
    // uses all-split mode which writes full and partial reads into separate
    // subdirectories per sample:
    //   <sample>/full/*.fastq.gz
    //   <sample>/partial/*.fastq.gz
    // For downstream processing (subsampling, alignment) we only want the
    // full reads. For publishing we want the entire sample directory so both
    // partial and full are available.
    // Duplex adapter types use full-only mode (flat structure), so the sample
    // directory is passed through unchanged. We detect the split layout by
    // checking for the full/ subdirectory rather than duplicating the adapter
    // type list from the config.
    ch_full_reads = ch_samples_with_meta.map { meta, sample_dir ->
        def full_dir = sample_dir.resolve('full')
        [meta, full_dir.exists() ? full_dir : sample_dir]
    }

    ch_demux_metrics_json_report = XOOS_PARSE_DEMUX_METRICS(
        demux_output.demux_metrics.map { meta, metrics_dir -> [meta, metrics_dir.resolve("run_stats.tsv"), metrics_dir.resolve("sample_stats.tsv")] }
    ).report

    ch_multiqc_files = XOOS_DEMUX_METRICS_TO_MULTIQC(
        demux_output.demux_metrics.map { meta, metrics_dir -> [meta, metrics_dir.resolve("run_stats.tsv"), metrics_dir.resolve("sample_stats.tsv")] }
    ).multiqc_report_dir

    // Collect versions (also published via topic channels)
    ch_versions = Channel.empty()
    ch_versions = ch_versions.mix(XOOS_DEMUX.out.versions)
    ch_versions = ch_versions.mix(XOOS_PARSE_DEMUX_METRICS.out.versions)
    ch_versions = ch_versions.mix(XOOS_DEMUX_METRICS_TO_MULTIQC.out.versions)

    emit:
    demultiplexed_fastqs     = ch_full_reads         // channel: [ [run_name: name, run_type: type, id: name, ...], path(full_reads_dir) ]
    demultiplexed_fastqs_all = ch_samples_with_meta  // channel: [ [run_name: name, run_type: type, id: name, ...], path(sample_dir) ]
    report                   = ch_demux_metrics_json_report // channel: [ val(meta), path(report) ]
    metrics_files            = demux_output.demux_metrics // channel: [ val(meta), path(metrics_dir) ]
    multiqc_files            = ch_multiqc_files // channel: [ val(meta), path(multiqc_report_dir) ]
    versions                 = ch_versions
}
