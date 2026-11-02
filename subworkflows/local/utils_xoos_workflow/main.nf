include { samplesheetToList } from 'plugin/nf-schema'
include { convertSamplesheet } from '../utils_samplesheet/main'

// ---------------------------------------------------------------------------
// R1: Helper functions extracted from workflows/xoosnf.nf
// ---------------------------------------------------------------------------

/*
 * Helper: remove a key from a map and return the new map
 */
def removeKey(map, key) {
    def new_map = map.clone()
    new_map.remove(key)
    return new_map
}

/*
 * Helper: map a pipeline run type to the library type expected by XOOS tools.
 */
def libraryTypeFromRunType(runType) {
    def normalizedRunType = runType?.toString()?.toLowerCase()
    def libraryType = [
        'sbx-d': 'sbx-d',
        'sbx-dm': 'sbx-d',
        'sbx-fast': 'sbx-fast',
        'i2x_hpsid_v0.01': 'sbx-d',
        'i2x_hpsid_v0.00': 'sbx-d',
        'i2x_hpu_v0.11': 'sbx-d',
        'i2x_hpu_v0.10': 'sbx-d',
        'ys': 'sbx-s',
        'ysu': 'sbx-s',
        'ys-new': 'sbx-s',
    ][normalizedRunType]

    if (libraryType == null) {
        throw new IllegalArgumentException("Unsupported run type '${runType}'. Expected SBX-D, SBX-DM, SBX-FAST, I2X_HPSID_v0.01, I2X_HPSID_v0.00, I2X_HPU_v0.11, I2X_HPU_v0.10, YS, YSU, or YS-NEW.")
    }
    return libraryType
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
    def samplesheet_file = convertSamplesheet(file(samplesheet_path, checkIfExists: true))
    samplesheetToList(samplesheet_file, "${projectDir}/assets/schema_samplesheet.json")
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
 * Resolve the static graph components from the filtered Giraffe index.
 * This is also used for pre-aligned BAM/CRAM inputs, where the graph used for
 * the original alignment is unknown and no per-sample graph tuple is emitted.
 */
def resolveGiraffeIndexComponents(idxDir) {
    def baseDir = file(idxDir, type: 'dir', checkIfExists: true)
    def gbzFiles = baseDir.list()?.findAll { it -> it.endsWith('.gbz') }?.sort()
    if (!gbzFiles) {
        error("No .gbz file found in Giraffe index directory: ${idxDir}")
    }
    def prefix = gbzFiles.first().replaceFirst(/\.gbz$/, '')

    return [
        gbz       : file(baseDir.resolve("${prefix}.gbz"), checkIfExists: true),
        zipcodes  : file(baseDir.resolve("${prefix}.zipcodes"), checkIfExists: true),
        minimizer : file(baseDir.resolve("${prefix}.min"), checkIfExists: true),
        dist      : file(baseDir.resolve("${prefix}.dist"), checkIfExists: true),
    ]
}

/*
 * Resolve the pre-generated chrM path GAF bundled with a Giraffe index.
 * The file uses the graph prefix followed by .chrm_paths.grch38.gaf.
 */
def resolveChrmPathsFromIndexDir(idxDir) {
    def baseDir = file(idxDir, type: 'dir', checkIfExists: true)
    def gbzFiles = baseDir.list()?.findAll { it -> it.endsWith('.gbz') }?.sort()
    if (!gbzFiles) {
        error("No .gbz file found in Giraffe index directory: ${idxDir}")
    }
    def prefix = gbzFiles.first().replaceFirst(/\.gbz$/, '')
    return file(baseDir.resolve("${prefix}.chrm_paths.grch38.gaf"), checkIfExists: true)
}

// ---------------------------------------------------------------------------
// R2: channel.value factory helpers — collapse 3-line resource channel
//     construction into single-line calls throughout the main workflow.
// ---------------------------------------------------------------------------

/*
 * Wrap an optional single file in a [meta, file] channel.value tuple.
 * Returns [[], []] when path is null/empty.
 */
def makeFileChannel(String id, path) {
    return path
        ? channel.value([[id: id], file(path, checkIfExists: true)])
        : channel.value([[], []])
}

/*
 * Wrap an optional directory in a [meta, dir] channel.value tuple.
 * Returns [[], []] when dirPath is null/empty.
 */
def makeDirChannel(String id, dirPath) {
    return dirPath
        ? channel.value([[id: id], file(dirPath, type: 'dir', checkIfExists: true)])
        : channel.value([[], []])
}

/*
 * Wrap optional paired files in a [meta, file1, file2] channel.value tuple.
 * Returns [[], [], []] when either path is null/empty.
 */
def makeFilePairChannel(String id, path1, path2) {
    return (path1 && path2)
        ? channel.value([[id: id], file(path1, checkIfExists: true), file(path2, checkIfExists: true)])
        : channel.value([[], [], []])
}

/*
 * Wrap an optional FASTA+FAI pair in a [meta, fasta, fai] channel.value tuple.
 * Returns [[], [], []] when fastaPath is null/empty.
 */
def makeFastaChannel(String id, fastaPath, faiPath) {
    return fastaPath
        ? channel.value([[id: id], file(fastaPath, checkIfExists: true), file(faiPath, checkIfExists: true)])
        : channel.value([[], [], []])
}

// ---------------------------------------------------------------------------
// R3: Variant-calling reference resolution
// ---------------------------------------------------------------------------

/*
 * Resolve variant-calling reference channels.
 *
 * Giraffe aligners map against the pangenome linearized reference rather than
 * the standard reference, so the VC reference must match. Also probes for an
 * optional <prefix>.subset.fa alongside the linearized reference which is used
 * as the partitioning reference (GATK4_SPLITINTERVALS) so that partitioning
 * only scans the subset chromosomes. Falls back to the full reference when no
 * subset file is present.
 *
 * Returns a map with keys: reference, referenceDict,
 * splitintervalsReference, splitintervalsReferenceDict.
 */
def resolveVcReferences(params, ch_reference) {
    def candidate_linearized_vc_ref = params.enable_pangenome_personalization.toBoolean()
        ? params.unfiltered_giraffe_reference_linearized
        : params.giraffe_reference_linearized
    def use_linearized_vc_ref = params.aligner in ['pb_giraffe', 'vg_giraffe'] && candidate_linearized_vc_ref

    def vc_reference = use_linearized_vc_ref
        ? channel.value([[id: 'vc_reference'], file(candidate_linearized_vc_ref, checkIfExists: true), file(appendExtension(candidate_linearized_vc_ref, '.fai'), checkIfExists: true)])
        : ch_reference
    def vc_reference_dict = use_linearized_vc_ref
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
    // variant callers themselves still use the full vc_reference. When the aligner is
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
    def splitintervals_reference = use_subset_ref
        ? channel.value([[id: 'splitintervals_reference'], candidate_subset_ref, file(appendExtension(candidate_subset_ref, '.fai'), checkIfExists: true)])
        : vc_reference
    def splitintervals_reference_dict = use_subset_ref
        ? channel.value([[id: 'splitintervals_reference_dict'], file(replaceExtension(candidate_subset_ref, '.dict'), checkIfExists: true)])
        : vc_reference_dict

    return [
        reference:                  vc_reference,
        referenceDict:              vc_reference_dict,
        splitintervalsReference:    splitintervals_reference,
        splitintervalsReferenceDict: splitintervals_reference_dict,
    ]
}

// ---------------------------------------------------------------------------
// R4: Somatic tumor-normal BAM pairing
// ---------------------------------------------------------------------------

/*
 * Group BAMs into tumor-normal pairs by the `group` samplesheet column.
 *
 * Emits one [meta, [tumor_bam, normal_bam], [tumor_bai, normal_bai]] tuple per
 * paired group. Singleton groups (unpaired BAMs) are logged as warnings and
 * dropped. Errors if no sample carries TN+group (which would indicate a
 * misconfigured somatic run).
 *
 * The pairing key is [group:] only — NOT the full meta — so tumor and normal
 * BAMs are grouped together correctly.
 */
def pairTumorNormalBams(ch_bams) {
    return ch_bams
        .filter { meta, _bam, _bai -> meta.containsKey('TN') && meta.containsKey('group') && meta.TN && meta.group }
        // error() is intentional here — do not downgrade to log.warn().
        // ifEmpty injects the closure's return value into the channel when empty;
        // error() throws, preventing a null from reaching the downstream .map/.groupTuple.
        .ifEmpty { error("variant_calling_mode is 'somatic_tn' but no samples have 'TN' and 'group' columns set in the samplesheet. Somatic tumor-normal variant calling requires paired samples with TN='T'/'N' and a shared group identifier.") }
        // The pairing key is deliberately reduced to [group:] only — it must
        // NOT use the full meta, which differs between tumor and normal samples.
        .map { meta, bam, bai -> [[group: meta.group], [meta, bam, bai]] }
        // remainder: true so singleton groups are emitted and can be logged below
        // rather than silently discarded.
        .groupTuple(size: 2, remainder: true)
        .filter { group, samples ->
            if (samples.size() < 2) {
                log.warn("Skipping group '${group.group}': expected a tumor-normal pair but found only ${samples.size()} sample(s). Check that both TN='T' and TN='N' are set for every sample in this group.")
                return false
            }
            return true
        }
        .map { group, samples ->
            def tumor = samples.find { meta, _bam, _bai -> meta.TN && meta.TN.toString().toUpperCase() == 'T' }
            def normal = samples.find { meta, _bam, _bai -> meta.TN && meta.TN.toString().toUpperCase() == 'N' }
            if (tumor && normal) {
                def meta = tumor[0].clone()
                meta.remove('TN')
                meta.remove('sample_name')
                meta.remove('sample_sid')
                meta.remove('group')
                meta.remove('bam_dedup_type')
                meta.remove('bam_secondary_rescued')
                meta.id = group.group
                meta.tumor_id = tumor[0].id
                meta.normal_id = normal[0].id

                // Missing or conflicting values use SVC's automatic sample-type detection.
                def tumor_sample_type = tumor[0].sample_type ?: 'auto-detect'
                def normal_sample_type = normal[0].sample_type ?: 'auto-detect'
                if (tumor_sample_type != normal_sample_type) {
                    log.warn("sample_type mismatch for group ${group.group}: tumor='${tumor_sample_type}', normal='${normal_sample_type}'. Defaulting to 'auto-detect'.")
                    meta.sample_type = 'auto-detect'
                } else {
                    meta.sample_type = tumor_sample_type
                }
                meta.tumor_run_name = tumor[0].run_name
                meta.normal_run_name = normal[0].run_name

                // Return [meta, [tumor_bam, normal_bam], [tumor_bai, normal_bai]]
                return [meta, [tumor[1], normal[1]], [tumor[2], normal[2]]]
            }
            else {
                log.warn("Could not identify tumor and normal BAM files for group ${group.group}.")
                null
            }
        }
        .filter { it -> it != null }
}


/*
 * Helper: resolve alignment_metrics_additional_regions labels against the
 * alignment_metrics_additional_region_beds map. Labels not found in the map are
 * skipped with a log.error. Returns a list of maps:
 * [ [region_name: '<label>', bed: '<path>'], ... ]
 */
def resolveAdditionalAlignmentMetricsRegions(regions_list, region_beds_map) {
    def resolved = []
    // Fail if the labels are not unique either in the regions_list or in the region_beds_map.
    if (regions_list.size() != regions_list.toSet().size()) {
        error "Duplicate labels found in regions_list: ${regions_list}"
    }
    if (region_beds_map.size() != region_beds_map.keySet().size()) {
        error "Duplicate labels found in region_beds_map: ${region_beds_map.keySet()}"
    }
    regions_list.each { label ->
        if (region_beds_map.containsKey(label)) {
            resolved << [region_name: label.toString(), bed: region_beds_map[label].toString()]
        } else {
            log.warn("alignment_metrics_additional_regions label '${label}' not found in alignment_metrics_additional_region_beds. Skipping.")
        }
    }
    return resolved
}


// Placeholder workflow required for Nextflow DSL2 module files
workflow UTILS_XOOS_WORKFLOW {
    main:
    versions = Channel.empty()
    emit:
    versions
}
