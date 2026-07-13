include { PREPARE_GIRAFFE_INPUTS } from './prepare_giraffe_inputs'
include { PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS } from '../../../modules/local/parabricks/giraffe/main'

workflow ALIGN_PB_GIRAFFE {
    take:
    ch_fastq // channel: [ val(meta), path(reads) ]
    ch_giraffe_index // channel: [ [id:'giraffe_idx'], path(idx_dir) ]
    ch_ref_paths // channel: [ path(dict) ] — linearized reference .dict for BAM @SQ header ordering

    main:
    giraffe_inputs = PREPARE_GIRAFFE_INPUTS(
        ch_fastq,
        ch_giraffe_index,
        params.enable_pangenome_personalization.toBoolean(),
    ).align_inputs

    pb_out = PARABRICKS_GIRAFFE_SINGLE_END_MULTIPLE_FASTQS(
        giraffe_inputs,
        ch_ref_paths,
    )
    bams = pb_out.bam

    emit:
    bams     = bams // channel: [ val(meta), path(bam), path(bai) ]
    versions = Channel.empty() // versions are collected via topic channels
}
