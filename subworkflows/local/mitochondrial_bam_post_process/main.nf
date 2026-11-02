include { XOOS_BAM_POST_PROCESS_RECRUIT        } from '../../../modules/local/xoos_bam_post_process/recruit/main'
include { VG_GIRAFFE_GAF                       } from '../../../modules/local/vg/giraffe_gaf/main'
include { XOOS_BAM_POST_PROCESS_CALL           } from '../../../modules/local/xoos_bam_post_process/call/main'
include { VG_SURJECT                                                                      } from '../../../modules/local/vg/surject/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_ORIGINAL_HEADER } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_ADDREPLACERG as SAMTOOLS_ADDREPLACERG_NFCORE } from '../../../modules/nf-core/samtools/addreplacerg/main'
include { SAMTOOLS_SORT_CHRM                                                                   } from './samtools_sort_chrm'

def readReadGroup(meta, header_file) {
    def read_group = header_file.text.readLines().find { it.startsWith('@RG') }
    if (!read_group) {
        if (workflow.stubRun) {
            return "'@RG\\tID:${meta.id}\\tSM:${meta.id}'"
        }
        error("No @RG header line found in ${header_file}.")
    }
    return "'${read_group.replace('\t', '\\t')}'"
}

workflow MITOCHONDRIAL_BAM_POST_PROCESS {
    take:
    ch_bams           // channel: [ val(meta), path(bam), path(bai) ]
    ch_giraffe_graphs // channel: [ val(meta), path(gbz), path(zipcodes), path(min), path(dist) ]
    ch_giraffe_chrm_paths // channel: path(pre-generated chrM paths GAF)

    main:
    ch_versions = Channel.empty()

    graph_components_by_sample = ch_giraffe_graphs
        .map { meta, graph, zipcodes, minimizers, dist ->
            [meta.id, graph, zipcodes, minimizers, dist]
        }

    XOOS_BAM_POST_PROCESS_RECRUIT(ch_bams)

    giraffe_inputs = XOOS_BAM_POST_PROCESS_RECRUIT.out.recruited
        .map { meta, _names, _regions, reads_fastq -> [meta.id, meta, reads_fastq] }
        .join(graph_components_by_sample, failOnMismatch: true, failOnDuplicate: true)
        .map { _sample_id, meta, reads_fastq, graph, zipcodes, minimizers, dist ->
            [meta, reads_fastq, graph, zipcodes, minimizers, dist]
        }

    VG_GIRAFFE_GAF(giraffe_inputs)

    call_inputs = ch_bams
        .join(XOOS_BAM_POST_PROCESS_RECRUIT.out.recruited)
        .join(VG_GIRAFFE_GAF.out.gaf)
        .map { meta, orig_bam, orig_bam_index, _names, recruit_regions, _reads_fastq, all_gaf ->
            [meta, all_gaf, orig_bam, orig_bam_index, recruit_regions]
        }
        .combine(ch_giraffe_chrm_paths)
        .map { meta, all_gaf, orig_bam, orig_bam_index, recruit_regions, chrm_paths_gaf ->
            [meta, all_gaf, chrm_paths_gaf, orig_bam, orig_bam_index, recruit_regions]
        }

    XOOS_BAM_POST_PROCESS_CALL(call_inputs)

    surject_inputs = ch_bams
        .join(XOOS_BAM_POST_PROCESS_CALL.out.gaf, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, orig_bam, _orig_bam_index, chrm_gaf ->
            [meta.id, meta, chrm_gaf]
        }
        .join(graph_components_by_sample, failOnMismatch: true, failOnDuplicate: true)
        .map { _sample_id, meta, chrm_gaf, graph, _zipcodes, _minimizers, _dist ->
            [meta, chrm_gaf, graph]
        }

    VG_SURJECT(surject_inputs, params.giraffe_mitochondrial_bam_post_process_chrm_name)

    add_rg_inputs = VG_SURJECT.out.bam
        .join(ch_bams, failOnMismatch: true, failOnDuplicate: true)
        .map { meta, unsorted_bam, orig_bam, orig_bam_index ->
            [meta, unsorted_bam, orig_bam, orig_bam_index]
        }

    ch_original_bams = add_rg_inputs.map { meta, _unsorted_bam, orig_bam, orig_bam_index ->
        [meta, orig_bam, orig_bam_index]
    }

    SAMTOOLS_VIEW_ORIGINAL_HEADER(
        ch_original_bams,
        [[], [], []],
        [[], []],
        [[], []],
        '',
    )

    ch_read_groups = SAMTOOLS_VIEW_ORIGINAL_HEADER.out.sam
        .map { meta, header_file ->
            [meta.id, readReadGroup(meta, header_file)]
        }

    ch_addreplacerg_inputs = add_rg_inputs
        .map { meta, unsorted_bam, _orig_bam, orig_bam_index ->
            [meta.id, meta, unsorted_bam, orig_bam_index]
        }
        .join(ch_read_groups, by: 0)
        .map { _sample_id, meta, unsorted_bam, orig_bam_index, read_group ->
            [meta, unsorted_bam, orig_bam_index, read_group]
        }

    SAMTOOLS_ADDREPLACERG_NFCORE(
        ch_addreplacerg_inputs,
        [[], [], [], []],
    )
    SAMTOOLS_SORT_CHRM(SAMTOOLS_ADDREPLACERG_NFCORE.out.bam)

    emit:
    bam      = SAMTOOLS_SORT_CHRM.out.bam // channel: [ val(meta), path(chrM_bam), path(chrM_bai) ]
    versions = ch_versions
}
