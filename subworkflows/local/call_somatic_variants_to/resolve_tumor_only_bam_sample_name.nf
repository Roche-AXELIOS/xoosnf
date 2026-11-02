include { SAMTOOLS_BAM_SAMPLE_NAME } from '../../../modules/local/samtools/bam_sample_name/main'

def readSampleName(sample_name_file) {
    def sample_name = sample_name_file.text.trim()
    if (!sample_name) {
        error("The BAM sample-name resolver produced an empty SM tag for ${sample_name_file}.")
    }
    return sample_name
}

workflow RESOLVE_TUMOR_ONLY_BAM_SAMPLE_NAME {
    take:
    ch_tumor_bams // channel: [meta, [bam_tumor], [bai_tumor]]

    main:
    ch_sample_name_inputs = ch_tumor_bams.map { meta, bams, _bais ->
        [meta.id, meta, 'tumor', bams[0]]
    }
    SAMTOOLS_BAM_SAMPLE_NAME(ch_sample_name_inputs)

    ch_sample_names = SAMTOOLS_BAM_SAMPLE_NAME.out.sample_name
        .map { sample_key, _meta, _role, sample_name_file ->
            [sample_key, readSampleName(sample_name_file)]
        }

    bams = ch_tumor_bams
        .map { meta, tumor_bams, tumor_bam_indices -> [meta.id, meta, tumor_bams, tumor_bam_indices] }
        .join(ch_sample_names, by: 0)
        .map { _sample_key, meta, tumor_bams, tumor_bam_indices, tumor_sample_name ->
            [meta + [tumor_sample_name: tumor_sample_name], tumor_bams, tumor_bam_indices]
        }

    emit:
    bams
}
