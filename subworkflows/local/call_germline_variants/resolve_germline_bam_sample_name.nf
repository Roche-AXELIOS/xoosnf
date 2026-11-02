include { SAMTOOLS_BAM_SAMPLE_NAME } from '../../../modules/local/samtools/bam_sample_name/main'

def readSampleName(sample_name_file) {
    def sample_name = sample_name_file.text.trim()
    if (!sample_name) {
        error("The BAM sample-name resolver produced an empty SM tag for ${sample_name_file}.")
    }
    return sample_name
}

workflow RESOLVE_GERMLINE_BAM_SAMPLE_NAME {
    take:
    ch_bams // channel: [meta, bam, bai]

    main:
    ch_sample_name_inputs = ch_bams.map { meta, bam, _bai ->
        [meta.id, meta, 'sample', bam]
    }
    SAMTOOLS_BAM_SAMPLE_NAME(ch_sample_name_inputs)

    ch_sample_names = SAMTOOLS_BAM_SAMPLE_NAME.out.sample_name
        .map { sample_key, _meta, _role, sample_name_file ->
            [sample_key, readSampleName(sample_name_file)]
        }

    bams = ch_bams
        .map { meta, bam, bai -> [meta.id, meta, bam, bai] }
        .join(ch_sample_names, by: 0)
        .map { _sample_key, meta, bam, bai, bam_sample_name ->
            [meta + [bam_sample_name: bam_sample_name], bam, bai]
        }

    emit:
    bams
}
