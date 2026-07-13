include { PREPARE_GIRAFFE_INPUTS } from './prepare_giraffe_inputs'
include { VG_GIRAFFE             } from '../../../modules/local/vg/giraffe/main'

workflow ALIGN_VG_GIRAFFE {
    take:
    ch_fastq // channel: [ val(meta), path(reads) ]
    ch_giraffe_index // channel: [ [id:'giraffe_idx'], path(idx_dir) ]

    main:
    giraffe_inputs = PREPARE_GIRAFFE_INPUTS(
        ch_fastq,
        ch_giraffe_index,
        params.enable_pangenome_personalization.toBoolean(),
    ).align_inputs

    bams = VG_GIRAFFE(
        giraffe_inputs
    ).bam

    emit:
    bams     = bams // channel: [ val(meta), path(bam), path(bai) ]
    versions = Channel.empty() // versions are collected via topic channels
}
