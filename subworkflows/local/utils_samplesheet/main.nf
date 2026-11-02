//
// Samplesheet format detection and conversion utilities.
//
// Supports two input formats:
// 1. Flat CSV — standard format with sample_name, sample_sid columns.
// 2. Sectioned — INI-style with a [Samples] header, Sample ID / Sample Index columns,
//    and optional other sections ([Header], [Settings], etc.) that are ignored.
//
// Comment lines (starting with #) and blank lines are stripped from both formats.
//

// The def functions below (convertSamplesheet, validateSamplesheet) are imported
// directly by pipeline workflows; UTILS_SAMPLESHEET itself is not called as a workflow.
//

/*
 * Check whether a line is a section header (e.g. [Samples], [Header]).
 * Only matches lines whose entire content is a bracketed name — lines containing
 * commas are excluded to avoid false positives on CSV data values like "[control]".
 */
def isSectionHeader(String trimmedLine) {
    return trimmedLine ==~ /(?i)\[[^\]]+\]/ && !trimmedLine.contains(',')
}

/*
 * Column name aliases: maps normalized names (lowercase, underscores) to canonical
 * schema property names. Includes all schema properties to preserve correct casing
 * (e.g., 'TN' must stay uppercase since the schema uses additionalProperties: false).
 */
def samplesheetColumnAliases() {
    return [
        'sample_id'       : 'sample_name',
        'sample_index'    : 'sample_sid',
        'sample_name'     : 'sample_name',
        'sample_sid'      : 'sample_sid',
        'truth_sample_name': 'truth_sample_name',
        'target_bed'      : 'target_bed',
        'capture_bed'     : 'capture_bed',
        'target_coverage' : 'target_coverage',
        'group'           : 'group',
        'tn'              : 'TN',
        'sample_type'     : 'sample_type',
        'reporter_vcf'    : 'reporter_vcf',
        'reporter_vcf_tumor_sample' : 'reporter_vcf_tumor_sample',
        'reporter_vcf_normal_sample': 'reporter_vcf_normal_sample',
        'bam_type'        : 'bam_type',
    ]
}

/*
 * Resolve the sample identifier used for truth resources and vcfeval.
 *
 * An explicit samplesheet value supports truth samples whose name cannot be
 * inferred from the library or tumor-normal group. Callers may supply their
 * existing workflow-specific fallback; otherwise the first underscore-delimited
 * token is used.
 */
def resolveTruthSampleName(meta, fallback_name = meta.id.toString().split('_')[0]) {
    def configured_name = meta.truth_sample_name?.toString()?.trim()
    return configured_name ?: fallback_name
}

workflow UTILS_SAMPLESHEET {
    main:
    versions = Channel.empty()

    emit:
    versions
}

/*
 * Validate param-dependent samplesheet constraints.
 *
 * Called after format conversion to catch configuration errors before any
 * processing begins. Checks:
 * - enable_subsample_bam requires target_coverage per sample (or a global fallback).
 * - variant_calling_mode 'somatic_tn' requires TN and group columns with values.
 */
def validateSamplesheet(samplesheet_path) {
    def rows = file(samplesheet_path).splitCsv(header: true)

    if (params.enable_subsample_bam.toBoolean() && params.subsample_bam_coverage_target == null) {
        rows.each { row ->
            if (!row.containsKey('target_coverage') || row.target_coverage == null || row.target_coverage.toString().trim() == '') {
                error "Sample '${row.sample_name}' has no target_coverage in the samplesheet. " +
                      "When enable_subsample_bam is true, either provide a 'target_coverage' column " +
                      "in the samplesheet or set params.subsample_bam_coverage_target as a global fallback."
            }
        }
    }

    if (params.variant_calling_mode == 'somatic_tn') {
        def hasAnyPair = rows.any { row ->
            row.containsKey('TN') && row.containsKey('group') &&
            row.TN != null && row.TN.toString().trim() != '' &&
            row.group != null && row.group.toString().trim() != ''
        }
        if (!hasAnyPair) {
            error "variant_calling_mode is 'somatic_tn' but no samples have 'TN' and 'group' columns " +
                  "set in the samplesheet. Somatic tumor-normal variant calling requires paired samples " +
                  "with TN='T'/'N' and a shared group identifier."
        }
    }
}



/*
 * Normalize a samplesheet for downstream consumption.
 *
 * Handles two concerns:
 * 1. Comment stripping — lines starting with # are removed from any samplesheet.
 * 2. Sectioned format conversion — files with a [Samples] section header are
 *    converted to canonical flat CSV (column aliases applied, other sections dropped).
 *
 * Returns the original path unchanged when no transformation is needed.
 */
def convertSamplesheet(samplesheet_path) {
    def samplesheet_file = file(samplesheet_path)
    def rawContent = samplesheet_file.text
    def lines = rawContent.readLines()

    def columnAliases = samplesheetColumnAliases()

    // Determine format: sectioned (has section headers like [Samples]) or flat CSV
    def hasSections = lines.any { line -> isSectionHeader(line.trim()) }

    // Check for comment lines or blank lines that need stripping
    def needsCleaning = lines.any { line ->
        def trimmed = line.trim()
        trimmed.startsWith('#') || trimmed == ''
    }

    if (hasSections) {
        return convertSectionedSamplesheet(samplesheet_file, rawContent, lines, columnAliases)
    }

    if (needsCleaning) {
        return stripCommentsAndBlanks(rawContent, lines)
    }

    // No transformation needed
    return samplesheet_file
}

/*
 * Strip comment lines and blank lines from a flat CSV samplesheet.
 */
def stripCommentsAndBlanks(String rawContent, lines) {
    def cleanedLines = lines.findAll { line ->
        def trimmed = line.trim()
        !trimmed.startsWith('#') && trimmed != ''
    }
    return writeTempSamplesheet(rawContent, cleanedLines)
}

/*
 * Convert a sectioned samplesheet to canonical flat CSV.
 * Extracts the [Samples] section, normalizes column names, applies aliases.
 */
def convertSectionedSamplesheet(samplesheet_file, String rawContent, lines, columnAliases) {
    // Find the [Samples] section
    def inSamplesSection = false
    def samplesFound = false
    def samplesLines = []
    lines.each { line ->
        def trimmed = line.trim()
        // Skip empty lines and comments
        if (trimmed == '' || trimmed.startsWith('#')) {
            return
        }
        // Check for section headers
        if (isSectionHeader(trimmed)) {
            if (trimmed ==~ /(?i)\[samples\]/) {
                inSamplesSection = true
                samplesFound = true
            } else {
                inSamplesSection = false
            }
            return
        }
        if (inSamplesSection) {
            samplesLines << trimmed
        }
    }

    if (!samplesFound) {
        error("Samplesheet '${samplesheet_file}' contains section headers but no [Samples] section. " +
              "A [Samples] section with 'Sample ID' and 'Sample Index' columns is required.")
    }

    if (samplesLines.size() < 2) {
        error("Samplesheet '${samplesheet_file}' has a [Samples] section but no data rows.")
    }

    // Parse header and normalize column names
    def rawHeader = samplesLines[0].split(',').collect { it.trim() }
    def expectedColumnCount = rawHeader.size()
    def canonicalHeader = rawHeader.collect { col ->
        def normalized = col.toLowerCase().replaceAll(/\s+/, '_')
        columnAliases.getOrDefault(normalized, normalized)
    }

    // Detect duplicate canonical column names (e.g. both 'Sample ID' and 'sample_name')
    def seen = [] as Set
    canonicalHeader.each { col ->
        if (!seen.add(col)) {
            error("Samplesheet '${samplesheet_file}' [Samples] section has duplicate column '${col}' " +
                  "after alias normalization. Original columns: ${rawHeader.join(', ')}")
        }
    }

    // Validate required columns are present
    if (!canonicalHeader.contains('sample_name')) {
        error("Samplesheet '${samplesheet_file}' [Samples] section is missing a 'Sample ID' column. " +
              "Found columns: ${rawHeader.join(', ')}")
    }
    if (!canonicalHeader.contains('sample_sid')) {
        error("Samplesheet '${samplesheet_file}' [Samples] section is missing a 'Sample Index' column. " +
              "Found columns: ${rawHeader.join(', ')}")
    }

    // Reorder columns so sample_name and sample_sid come first (matching schema order).
    // The xoos_demux stub reads the first two CSV columns positionally.
    def requiredFirst = ['sample_name', 'sample_sid']
    def remaining = canonicalHeader.findAll { !requiredFirst.contains(it) }
    def orderedHeader = requiredFirst + remaining

    // Build index mapping from new column order to original positions
    def orderIndices = orderedHeader.collect { col -> canonicalHeader.indexOf(col) }

    // Build the output CSV, validating column counts and reordering values
    def outputLines = [orderedHeader.join(',')]
    (1..<samplesLines.size()).each { i ->
        def values = samplesLines[i].split(',', -1).collect { it.trim() }
        if (values.size() != expectedColumnCount) {
            error("Samplesheet '${samplesheet_file}' [Samples] data row ${i} has ${values.size()} " +
                  "columns but the header has ${expectedColumnCount}.")
        }
        def reordered = orderIndices.collect { idx -> values[idx] }
        outputLines << reordered.join(',')
    }

    return writeTempSamplesheet(rawContent, outputLines)
}

/*
 * Write samplesheet lines to a deterministic temp file.
 * Uses a hash of the original file content to avoid collisions and stale cache.
 * Skips the write when the target already exists (idempotent).
 * Uses atomic write-then-rename to prevent parallel runs from reading a partial file.
 */
def writeTempSamplesheet(String rawContent, outputLines) {
    def content = outputLines.join('\n') + '\n'
    def hash = rawContent.md5()
    def convertedDir = java.nio.file.Paths.get(System.getProperty('java.io.tmpdir'), 'xoos-converted-samplesheets')
    java.nio.file.Files.createDirectories(convertedDir)
    def convertedPath = convertedDir.resolve("${hash}.csv")
    if (!java.nio.file.Files.exists(convertedPath)) {
        def tmpPath = convertedDir.resolve("${hash}.csv.tmp.${Thread.currentThread().id}")
        tmpPath.toFile().text = content
        java.nio.file.Files.move(tmpPath, convertedPath,
            java.nio.file.StandardCopyOption.ATOMIC_MOVE,
            java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    }
    return file(convertedPath)
}
