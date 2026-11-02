include { SAMTOOLS_BAM_SAMPLE_NAME } from '../../../modules/local/samtools/bam_sample_name/main'

def readSampleName(sample_name_file) {
    def sample_name = sample_name_file.text.trim()
    if (!sample_name) {
        error("The BAM sample-name resolver produced an empty SM tag for ${sample_name_file}.")
    }
    return sample_name
}

workflow RESOLVE_TUMOR_NORMAL_BAM_SAMPLE_NAMES {
    take:
    ch_tumor_normal_bams // channel: [meta, [bam_tumor, bam_normal], [bai_tumor, bai_normal]]

    main:
    ch_sample_name_inputs = ch_tumor_normal_bams.flatMap { meta, bams, _bais ->
        def inputs = [[meta.id, meta, 'tumor', bams[0]]]
        if (bams.size() > 1) {
            inputs << [meta.id, meta, 'normal', bams[1]]
        }
        inputs
    }
    SAMTOOLS_BAM_SAMPLE_NAME(ch_sample_name_inputs)

    ch_sample_names = SAMTOOLS_BAM_SAMPLE_NAME.out.sample_name
        .map { sample_key, _meta, role, sample_name_file ->
            [sample_key, role, readSampleName(sample_name_file)]
        }
    ch_tumor_sample_names = ch_sample_names
        .filter { _sample_key, role, _sample_name -> role == 'tumor' }
        .map { sample_key, _role, sample_name -> [sample_key, sample_name] }
    ch_normal_sample_names = ch_sample_names
        .filter { _sample_key, role, _sample_name -> role == 'normal' }
        .map { sample_key, _role, sample_name -> [sample_key, sample_name] }

    bams = ch_tumor_normal_bams
        .map { meta, tumor_normal_bams, tumor_normal_bam_indices -> [meta.id, meta, tumor_normal_bams, tumor_normal_bam_indices] }
        .join(ch_tumor_sample_names, by: 0)
        .join(ch_normal_sample_names, by: 0)
        .map { _sample_key, meta, tumor_normal_bams, tumor_normal_bam_indices, tumor_sample_name, normal_sample_name ->
            [meta + [tumor_sample_name: tumor_sample_name, normal_sample_name: normal_sample_name], tumor_normal_bams, tumor_normal_bam_indices]
        }

    emit:
    bams
}
