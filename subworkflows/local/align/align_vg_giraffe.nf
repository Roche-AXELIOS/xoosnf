include { VG_GIRAFFE             } from '../../../modules/local/vg/giraffe/main'

workflow ALIGN_VG_GIRAFFE {
    take:
    ch_giraffe_inputs // channel: [ val(meta), path(reads), path(graph), path(zipcodes), path(minimizers), path(dist) ]
    ch_ref_paths // channel: [ path(dict) ] — linearized reference .dict for vg giraffe --ref-paths

    main:
    vg_out = VG_GIRAFFE(
        ch_giraffe_inputs,
        ch_ref_paths,
    )
    bams = vg_out.bam

    emit:
    bams           = bams // channel: [ val(meta), path(bam), path(bai) ]
    command_args   = vg_out.command_args // channel: [ val(meta), val(args) ]
    versions       = Channel.empty() // versions are collected via topic channels
}
