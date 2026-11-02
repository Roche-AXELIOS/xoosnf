include { XOOS_SMALL_VARIANT_CALLER_COMPUTE_VCF_FEATURES         } from '../../../modules/local/xoos_small_variant_caller/compute_vcf_features/main'
include { XOOS_SMALL_VARIANT_CALLER_COMPUTE_BAM_FEATURES         } from '../../../modules/local/xoos_small_variant_caller/compute_bam_features/main'
include { XOOS_SMALL_VARIANT_CALLER_TRAIN_MODEL                  } from '../../../modules/local/xoos_small_variant_caller/train_model/main'
include { XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS } from '../../../modules/local/xoos_small_variant_caller/filter_variants_germline_wgs/main'

/*
 * Small variant caller with sample-specific model retraining
 * Orchestrates: compute features → train model → filter variants
 */
workflow RETRAIN_GERMLINE_SMALL_VARIANT_MODEL {
    take:
    ch_input // channel: [ val(meta), path(bams), path(bam_indexes), path(vcf), path(vcf_index) ]
    ch_reference // value: [ path(reference_fasta), path(reference_fasta_index) ]
    ch_par_beds // channel: [ path(chrX_par_bed), path(chrY_par_bed) ] - parameter BED for filtering (e.g. high-confidence regions)
    ch_popaf // channel: [ path(vcf.gz), path(vcf.gz.tbi) ]
    ch_train_bed // channel: [ val(meta), path(bed) ] - training regions (per sample)
    ch_test_bed // channel: [ val(meta), path(bed) ] - evaluation regions (per sample)
    ch_truth_vcf // channel: [ val(meta), path(vcf.gz), path(vcf.gz.tbi) ]
    aligner // value: aligner name (e.g. 'pb_giraffe', 'bwa')

    main:

    // Step 1: Compute VCF features (extract features from input VCF on training regions)
    XOOS_SMALL_VARIANT_CALLER_COMPUTE_VCF_FEATURES(
        ch_input.map { meta, _bams, _bam_indexes, vcf, vcf_index -> [meta, vcf, vcf_index] }.join(ch_train_bed),
        ch_reference,
        ch_popaf,
    )

    // Step 2: Compute BAM features (extract features from BAMs on VCF regions)
    XOOS_SMALL_VARIANT_CALLER_COMPUTE_BAM_FEATURES(
        ch_input.map { meta, bams, bam_indexes, _vcf, _vcf_index -> [meta, bams, bam_indexes] }.join(
            XOOS_SMALL_VARIANT_CALLER_COMPUTE_VCF_FEATURES.out.vcf_bed
        ),
        ch_reference,
    )

    // Step 3: Train models using computed features and truth VCF
    XOOS_SMALL_VARIANT_CALLER_TRAIN_MODEL(
        XOOS_SMALL_VARIANT_CALLER_COMPUTE_BAM_FEATURES.out.bam_features.join(
            XOOS_SMALL_VARIANT_CALLER_COMPUTE_VCF_FEATURES.out.vcf_features
        ).join(ch_truth_vcf)
    )

    // Step 4: Filter variants using trained models (on eval regions)
    // Join the per-sample inputs, models, and test BED by meta so each task
    // receives the resources for its own sample. These are independent queue
    // channels emitted in non-deterministic order by upstream processes; passing
    // them as separate process inputs would zip them positionally and mis-pair
    // samples (e.g. HG002 calls filtered against HG005 regions). Join on the
    // original meta first, then apply the downstream metadata in lockstep.
    ch_models = XOOS_SMALL_VARIANT_CALLER_TRAIN_MODEL.out.model_snv.join(
        XOOS_SMALL_VARIANT_CALLER_TRAIN_MODEL.out.model_indel)

    ch_filter_joined = ch_input
        .join(ch_models, failOnMismatch: true, failOnDuplicate: true)
        .join(ch_test_bed, failOnMismatch: true, failOnDuplicate: true)
        .multiMap { meta, bams, bam_indexes, vcf, vcf_index, model_snv, model_indel, test_bed ->
            def filter_meta = meta + [model_type: 'retrained', vcfeval_bed: 'giabhc_chr4-22']
            input: [filter_meta, bams, bam_indexes, vcf, vcf_index]
            models: [filter_meta, model_snv, model_indel]
            test_bed: [filter_meta, test_bed]
        }

    XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS(
        ch_filter_joined.input,
        ch_reference,
        ch_filter_joined.models,
        ch_popaf,
        ch_par_beds,
        ch_filter_joined.test_bed,
        aligner,
    )

    emit:
    vcf         = XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS.out.vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    metrics     = XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS.out.metrics // channel: [ val(meta), path(metrics_tsv) ]
    models      = XOOS_SMALL_VARIANT_CALLER_TRAIN_MODEL.out.model_snv.join(XOOS_SMALL_VARIANT_CALLER_TRAIN_MODEL.out.model_indel) // channel: [ val(meta), path(snv), path(indel) ]
}
