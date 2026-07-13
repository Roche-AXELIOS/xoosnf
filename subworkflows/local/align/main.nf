// Import aligner-specific subworkflows
include { ALIGN_BWA          } from './align_bwa'
include { ALIGN_MINIMAP2     } from './align_minimap2'
include { ALIGN_PB_FQ2BAM    } from './align_pb_fq2bam'
include { ALIGN_PB_MINIMAP2  } from './align_pb_minimap2'
include { ALIGN_PB_GIRAFFE   } from './align_pb_giraffe'
include { ALIGN_VG_GIRAFFE   } from './align_vg_giraffe'

workflow ALIGN {
    take:
    aligner          // val: 'bwa', 'minimap2', etc.
    ch_fastq         // channel: [ val(meta), path(reads) ]
    ch_reference     // channel: [ [id:'reference'], path(fasta) ]
    ch_bwa_idx       // channel: [ [id:'bwa_idx'], path(idx_dir) ]
    ch_giraffe_idx   // channel: [ [id:'giraffe_idx'], path(idx_dir) ]
    ch_ref_paths     // channel: [ path(dict) ] — linearized reference .dict for pbrun giraffe --ref-paths

    main:
    // Collect versions (also published via topic channels)
    ch_versions = Channel.empty()

    // 1. Handle conditional subsampling
    ch_reads_for_alignment = ch_fastq

    // 2. Initialize empty output channel
    bams = channel.empty()

    // 3. Select aligner and delegate to appropriate subworkflow using the conditionally modified channel
    if (aligner == 'bwa') {
        bwa_out = ALIGN_BWA(
            ch_reads_for_alignment,
            ch_reference,
            ch_bwa_idx,
        )
        bams = bwa_out.bams
        ch_versions = ch_versions.mix(ALIGN_BWA.out.versions)
    }
    else if (aligner == 'minimap2') {
        mm_out = ALIGN_MINIMAP2(
            ch_reads_for_alignment,
            ch_reference,
        )
        bams = mm_out.bams
        ch_versions = ch_versions.mix(ALIGN_MINIMAP2.out.versions)
    }
    else if (aligner == 'pb_fq2bam') {
        pb_out = ALIGN_PB_FQ2BAM(
            ch_reads_for_alignment,
            ch_reference,
            ch_bwa_idx,
        )
        bams = pb_out.bams
        ch_versions = ch_versions.mix(ALIGN_PB_FQ2BAM.out.versions)
    }
    else if (aligner == 'pb_minimap2') {
        pb_mm_out = ALIGN_PB_MINIMAP2(
            ch_reads_for_alignment,
            ch_reference,
        )
        bams = pb_mm_out.bams
        ch_versions = ch_versions.mix(ALIGN_PB_MINIMAP2.out.versions)
    }
    else if (aligner == 'pb_giraffe') {
        giraffe_out = ALIGN_PB_GIRAFFE(
            ch_reads_for_alignment,
            ch_giraffe_idx,
            ch_ref_paths,
        )
        bams = giraffe_out.bams
        ch_versions = ch_versions.mix(ALIGN_PB_GIRAFFE.out.versions)
    }
    else if (aligner == 'vg_giraffe') {
        giraffe_out = ALIGN_VG_GIRAFFE(
            ch_reads_for_alignment,
            ch_giraffe_idx,
        )
        bams = giraffe_out.bams
        ch_versions = ch_versions.mix(ALIGN_VG_GIRAFFE.out.versions)
    }
    else {
        error("Unsupported aligner specified: ${aligner}")
    }

    emit:
    bams     // channel: [ val(meta), path(bam), path(bai) ]
    versions = ch_versions
}
