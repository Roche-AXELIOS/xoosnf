include { GATK4_MUTECT2 as GATK4_MUTECT2_MITOCHONDRIAL                                            } from '../../../modules/local/gatk4/mutect2/main'
include { PARABRICKS_MUTECTCALLER as PARABRICKS_MUTECTCALLER_MITOCHONDRIAL                        } from '../../../modules/local/parabricks/mutectcaller/main'
include { XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_MITOCHONDRIA as XOOS_SVC_FILTER_MITOCHONDRIAL } from '../../../modules/local/xoos_small_variant_caller/filter_variants_mitochondria/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_MITOCHONDRIAL                                    } from '../../../modules/nf-core/bcftools/reheader/main'
include { MITOCHONDRIAL_BAM_POST_PROCESS                                                          } from '../mitochondrial_bam_post_process/main'

workflow MITOCHONDRIAL_VARIANT_CALLING {
    take:
    ch_bams           // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_reference_dict // channel: [ val(meta), path(dict) ]
    use_gatk_mutect2_mitochondrial  // val: boolean — when true use GATK Mutect2, when false use Parabricks mutectcaller
    ch_giraffe_graphs // channel: [ val(meta), path(gbz), path(zipcodes), path(min), path(dist) ]
    ch_giraffe_chrm_paths // channel: path(pre-generated chrM paths GAF)

    main:
    ch_versions = Channel.empty()
    ch_mutect_out = Channel.empty()
    ch_variant_bams = ch_bams

    // Match the pipeline-level default while preserving an explicit false.
    def enable_giraffe_mitochondrial_bam_post_process = params.enable_giraffe_mitochondrial_bam_post_process != null
        ? params.enable_giraffe_mitochondrial_bam_post_process.toBoolean()
        : (params.enable_mitochondrial_variant_caller.toBoolean() && params.aligner in ['pb_giraffe', 'vg_giraffe'])

    if (enable_giraffe_mitochondrial_bam_post_process) {
        MITOCHONDRIAL_BAM_POST_PROCESS(
            ch_bams,
            ch_giraffe_graphs,
            ch_giraffe_chrm_paths,
        )
        ch_variant_bams = MITOCHONDRIAL_BAM_POST_PROCESS.out.bam
    }

    if (use_gatk_mutect2_mitochondrial) {
        // Step 1a: Run GATK Mutect2 in mitochondria-mode on the configured
        // mitochondrial contig. No interval partitioning is needed because
        // the mitochondrial contig is a single ~16 kb contig.
        // The interval is passed via ext.args in nextflow.config.
        // Mutect2 emits a realigned BAM (--bamout) required by SVC for feature extraction.
        ch_mutect2_input = ch_variant_bams
            .map { meta, bam, bai -> [meta, [bam], [bai], []] }

        GATK4_MUTECT2_MITOCHONDRIAL(
            ch_mutect2_input,
            ch_reference.map { _meta, fasta, _fai -> [[id: 'fasta'], fasta] },
            ch_reference.map { _meta, _fasta, fai -> [[id: 'fai'], fai, []] },
            ch_reference_dict.map { _meta, dict -> [[id: 'dict'], dict] },
            channel.value([]),  // alleles
            channel.value([]),  // alleles_tbi
            channel.value([]),  // germline_resource
            channel.value([]),  // germline_resource_tbi
            channel.value([]),  // panel_of_normals
            channel.value([]),  // panel_of_normals_tbi
        )

        // Normalize: [meta, bam, bai, vcf, tbi] — join declared bai emit to avoid string-derived paths
        ch_mutect_out = GATK4_MUTECT2_MITOCHONDRIAL.out.bam
            .join(GATK4_MUTECT2_MITOCHONDRIAL.out.bai)
            .join(GATK4_MUTECT2_MITOCHONDRIAL.out.vcf)
            .join(GATK4_MUTECT2_MITOCHONDRIAL.out.tbi)
            .map { meta, bam, bai, vcf, tbi ->
                [meta, [bam], [bai], vcf, tbi]
            }
    } else {
        // Step 1b: Run Parabricks mutectcaller in tumor-only mode on the
        // configured mitochondrial contig. The interval and mitochondria-mode
        // are passed via ext.args in nextflow.config.
        // We reuse the standard module with empty normal BAM inputs.
        ch_pb_input = ch_variant_bams
            .map { meta, bam, bai -> [meta + [tumor_id: meta.id], bam, bai, [], [], []] }

        PARABRICKS_MUTECTCALLER_MITOCHONDRIAL(
            ch_pb_input,
            ch_reference.map { _meta, fasta, _fai -> [[id: 'fasta'], fasta] },
            [],  // panel_of_normals
            [],  // panel_of_normals_index
        )

        // Normalize: [meta, [bam], [bai], vcf, tbi]
        ch_mutect_out = PARABRICKS_MUTECTCALLER_MITOCHONDRIAL.out.bam
            .join(PARABRICKS_MUTECTCALLER_MITOCHONDRIAL.out.vcf)
            .map { meta, bam, bai, vcf, tbi ->
                [meta - meta.subMap('tumor_id'), [bam], [bai], vcf, tbi]
            }
    }

    // Step 2: Filter chrM variants with SVC germline-mitochondria workflow.
    // Uses AF/DP threshold filtering (no ML model required).
    // Input: Mutect2 realigned BAM + Mutect2 VCF
    XOOS_SVC_FILTER_MITOCHONDRIAL(
        ch_mutect_out,
        ch_reference,
    )

    // Step 3: Reheader the SVC VCF to set the sample name to meta.id.
    ch_mito_sample_names = XOOS_SVC_FILTER_MITOCHONDRIAL.out.vcf.collectFile { meta, _vcf, _tbi ->
        ["${meta.id}.txt", "${meta.id}" + '\n']
    }
    ch_mito_reheader_input = XOOS_SVC_FILTER_MITOCHONDRIAL.out.vcf
        .combine(ch_mito_sample_names)
        .filter { meta, _vcf, _tbi, sample_file -> sample_file.baseName == meta.id }
        .map { meta, vcf, _tbi, sample_file -> [meta, vcf, [], sample_file] }
    BCFTOOLS_REHEADER_MITOCHONDRIAL(
        ch_mito_reheader_input,
        channel.value([[], []]),
    )
    ch_reheadered_vcf = BCFTOOLS_REHEADER_MITOCHONDRIAL.out.vcf.join(BCFTOOLS_REHEADER_MITOCHONDRIAL.out.index)

    emit:
    vcf      = ch_reheadered_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    metrics  = XOOS_SVC_FILTER_MITOCHONDRIAL.out.metrics // channel: [ val(meta), path(metrics_tsv) ]
    versions = ch_versions
}
