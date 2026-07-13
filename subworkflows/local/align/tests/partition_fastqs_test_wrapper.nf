/*
 * Thin wrapper workflow to make partitionFastqs testable via nf-test.
 * partitionFastqs is a channel-operator function, not a workflow, so it
 * cannot be tested directly with nextflow_workflow or nextflow_function.
 */

include { partitionFastqs } from '../partition_fastqs'

workflow PARTITION_FASTQS_TEST {
    take:
    ch_fastq         // channel: [ val(meta), path(reads) ]
    partition_count  // val: number of partitions

    main:
    ch_partitioned = partitionFastqs(ch_fastq, partition_count)

    emit:
    partitions = ch_partitioned // channel: [ val(meta), list(reads) ]
}
