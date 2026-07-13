include { XOOS_ALIGNMENT_METRICS_COVERAGE } from '../../../modules/local/xoos_alignment_metrics/coverage/main'
include { SAMTOOLS_VIEW                   } from '../../../modules/nf-core/samtools/view/main'

/*
 * Helper: Parse a CSV file while skipping comment lines (starting with '#') and using the first non-comment line as the header for keys.
 */
def splitCsvSkipHeaderAndComment(file, sep) {
    def lines = file.readLines()
    def data = []
    // Iterate through lines, store header to use as keys, and skip comment lines starting with '#'
    def header = null
    lines.each { line ->
        if (line.startsWith('#')) {
            // Skip comment lines
            return null
        }
        if (header == null) {
            // First non-comment line is the header
            header = line.split(sep)
        }
        else {
            // Split the line and create a map using the header as keys
            def values = line.split(sep)
            def row = [:]
            header.eachWithIndex { key, index ->
                row[key] = values[index]
            }
            data << row
        }
    }
    return data
}

workflow SUBSAMPLE_BAM {
    take:
    ch_bam_bai // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta // channel: [ val(meta), path(fasta), path(fai) ]
    ch_bed // channel: [ val(meta), path(bed) ]
    target_coverage // val: val(depth) - fallback if target_coverage is not provided in meta
    metrics_name // val: tuple(val(metrics_row_name), val(metrics_column_name))

    main:
    // Collect versions (also published via topic channels)
    ch_versions = Channel.empty()

    // Compute coverage metrics using XOOS alignment_metrics
    ch_coverage_metrics = XOOS_ALIGNMENT_METRICS_COVERAGE(ch_bam_bai, ch_fasta, ch_bed).coverage_metrics
    ch_versions = ch_versions.mix(XOOS_ALIGNMENT_METRICS_COVERAGE.out.versions)
    // Parse the coverage metrics at 'coverage/coverage_stats.tsv' and extract the coverage by metrics name
    ch_coverage = ch_coverage_metrics
        .map { meta, path -> tuple(meta, path.resolve("coverage_stats.tsv")) }
        .map { meta, coverage_stats ->
            def coverage = 0.0
            // Skip one line when parsing CSV to ignore the comment line
            splitCsvSkipHeaderAndComment(file(coverage_stats), '\t').each { row ->
                if (row['metric_name'] == metrics_name[0]) {
                    coverage = row[metrics_name[1]] ? row[metrics_name[1]] as Integer : 0.0
                }
            }
            return tuple(meta, coverage)
        }
    // Output: [ val(meta), val(coverage) ]
    ch_ratio = ch_coverage.map { meta, coverage ->
        def sample_target = meta.containsKey('target_coverage') ? meta.target_coverage : target_coverage
        if (sample_target == null) {
            error "Sample '${meta.id}' has no target_coverage: provide it in the samplesheet or set params.subsample_bam_coverage_target"
        }
        def ratio = (coverage == 0) ? 1.0 : (sample_target as Integer) / coverage as Float
        return tuple(meta, Math.min(ratio, 1.0))
    }
    // Output: [ val(meta), val(ratio) ]
    // Combine the ratio with the original bam and bai to pass to SAMTOOLS_VIEW for downsampling
    ch_input_subsample = ch_bam_bai
        .join(ch_ratio)
        .map { meta, bam, bai, ratio ->
            return tuple(meta + ['subsample_fraction': ratio], bam, bai)
        }
    SAMTOOLS_VIEW(
        ch_input_subsample,
        [[], [], []],
        [[], []],
        [[], []],
        'bai',
    )
    ch_bam_subsampled = SAMTOOLS_VIEW.out.bam
        .join(
            SAMTOOLS_VIEW.out.bai
        )
    emit:
    bam_subsampled = ch_bam_subsampled // channel: [ val(meta), path(bam), path(bai) ]
    versions       = ch_versions
}
