include { SUBSAMPLE_RUN as SUBSAMPLE_BASECALL_RDB   } from '../../../modules/local/subsample/run/main'
include { SUBSAMPLE_RUN as SUBSAMPLE_BASECALL_FASTQ } from '../../../modules/local/subsample/run/main'

workflow SUBSAMPLE_BASECALL_RUNS {
    take:
    ch_basecall_rdb   // channel: [ val(meta), path(rdb_dir), path(samplesheet) ]
    ch_basecall_fastq // channel: [ val(meta), path(fastq_dir), path(samplesheet) ]

    main:
    // Collect versions (also published via topic channels)
    ch_versions = Channel.empty()

    ch_final_rdb   = ch_basecall_rdb
    ch_final_fastq = ch_basecall_fastq

    if (params.enable_subsample_run.toBoolean()) {
        // Strip the samplesheet off, run the module, and attach it back

        // Execute the RDB aliased module and update the channel
        SUBSAMPLE_BASECALL_RDB(ch_basecall_rdb)
        ch_final_rdb = SUBSAMPLE_BASECALL_RDB.out.subsampled
        ch_versions = ch_versions.mix(SUBSAMPLE_BASECALL_RDB.out.versions)

        // Execute the FASTQ aliased module and update the channel
        SUBSAMPLE_BASECALL_FASTQ(ch_basecall_fastq)
        ch_final_fastq = SUBSAMPLE_BASECALL_FASTQ.out.subsampled
        ch_versions = ch_versions.mix(SUBSAMPLE_BASECALL_FASTQ.out.versions)
    }

    ch_merged_basecalls = ch_final_rdb.mix(ch_final_fastq)

    emit:
    basecalls = ch_merged_basecalls // channel: [ val(meta), path(dir), path(samplesheet) ]
    versions  = ch_versions
}
