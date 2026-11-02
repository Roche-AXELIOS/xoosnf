include { PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS } from '../../../modules/local/parabricks/giraffe/main'

workflow ALIGN_PB_GIRAFFE {
    take:
    ch_giraffe_inputs // channel: [ val(meta), path(reads), path(graph), path(zipcodes), path(minimizers), path(dist) ]
    ch_ref_paths // channel: [ path(dict) ] — linearized reference .dict for BAM @SQ header ordering

    main:
    pb_out = PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS(
        ch_giraffe_inputs,
        ch_ref_paths,
    )
    bams = pb_out.bam

    emit:
    bams           = bams // channel: [ val(meta), path(bam), path(bai) ]
    command_args   = pb_out.command_args // channel: [ val(meta), val(args) ]
    versions       = Channel.empty() // versions are collected via topic channels
}
