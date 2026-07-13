/*
 * Partition a sample's FASTQ files into N groups for parallel alignment.
 *
 * Files are distributed using size-balanced bin-packing: sorted by size
 * descending, each file is assigned to the partition with the smallest
 * current total. This keeps total bytes per partition roughly equal
 * without decompressing or splitting any individual file.
 *
 * Always emits exactly partition_count items per sample. If a sample has
 * fewer files than partitions, the extra partitions are empty (file list
 * is []) and a warning is logged. Downstream code must filter or handle
 * empty partitions before alignment.
 *
 * Each emitted meta has:
 *   - id:          "${original_id}_${partition_index}" (unique process tag)
 *   - sample_id:   the true sample name (for read groups)
 *   - original_id: the pre-partition id (restored after merge)
 *   - partition:    0-based partition index
 *
 * Handles two input shapes:
 *   - [ meta, [file1.fastq.gz, file2.fastq.gz, ...] ]  (pre-demultiplexed)
 *   - [ meta, path(sample_dir) ]                        (demux output directory)
 */

def partitionFastqs(ch_fastq, partition_count) {
    // Nextflow 26.04 may pass workflow val inputs as DataflowVariable
    // wrappers.  The range operator (0..<n) requires a plain Integer.
    def n_partitions = partition_count as int

    return ch_fastq.flatMap { meta, reads ->
        // Resolve the list of FASTQ files from either a directory or a file list
        def files = []
        if (reads instanceof List) {
            files = reads.collect()
        } else {
            // Use file() with type:'dir' so Nextflow resolves through the
            // correct filesystem provider (local, S3, GCS). Path.list()
            // delegates to that provider, unlike File.eachFileRecurse()
            // which only works on local filesystems.
            def p = file(reads, checkIfExists: true)
            if (p.isDirectory()) {
                files = p.list()
                    .findAll { it.endsWith('.fastq.gz') }
                    .collect { file(p.resolve(it)) }
            } else {
                files = [reads]
            }
        }

        if (files.size() < n_partitions) {
            log.warn("Sample '${meta.id}': only ${files.size()} FASTQ file(s) available but ${n_partitions} partitions requested. ${n_partitions - files.size()} partition(s) will be empty.")
        }

        // Cache file sizes to avoid redundant stat calls (important for remote storage)
        def file_sizes = files.collectEntries { f -> [f, f.size()] }

        // Sort files by size descending for greedy bin-packing (toSorted
        // returns a new list, avoiding in-place mutation of the source)
        def sorted_files = files.toSorted { a, b -> file_sizes[b] <=> file_sizes[a] }

        // Always create exactly n_partitions partitions (some may be empty)
        def partitions = (0..<n_partitions).collect { [0L, []] }

        // Greedy assignment: put each file into the partition with the smallest total
        sorted_files.each { f ->
            def smallest = partitions.min { it[0] }
            smallest[0] += file_sizes[f]
            smallest[1] << f
        }

        // Emit exactly partition_count items per sample (empty partitions included)
        // Store original id so it can be restored after merge, and preserve
        // sample_id for read group construction (meta.id gets a _N suffix).
        partitions.withIndex()
            .collect { partition, idx ->
                def new_meta = meta.clone()
                new_meta.sample_id = meta.id
                new_meta.original_id = meta.id
                new_meta.partition = idx
                new_meta.id = "${meta.id}_${idx}"
                [new_meta, partition[1]]
            }
    }
}
