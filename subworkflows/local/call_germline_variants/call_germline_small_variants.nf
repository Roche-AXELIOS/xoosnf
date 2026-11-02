include { HAPLOTYPECALLER_GATK                                   } from './haplotypecaller_gatk'
include { HAPLOTYPECALLER_PARABRICKS                             } from './haplotypecaller_parabricks'
include { RETRAIN_GERMLINE_SMALL_VARIANT_MODEL                   } from './retrain_germline_small_variant_caller'

include { XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS } from '../../../modules/local/xoos_small_variant_caller/filter_variants_germline_wgs/main'

workflow CALL_GERMLINE_SMALL_VARIANTS {
    take:
    ch_bams // channel: [ val(meta), path(bam), path(bai) ]
    ch_reference // channel: [ val(meta), path(fasta), path(fai) ]
    ch_reference_dict // channel: [ val(meta), path(dict) ]
    ch_splitintervals_reference // channel: [ val(meta), path(fasta), path(fai) ] (subset reference for partitioning)
    ch_splitintervals_reference_dict // channel: [ val(meta), path(dict) ] (subset reference dict for partitioning)
    ch_par_beds // channel: [ val(meta), path(chrX_par_bed), path(chrY_par_bed) ]
    ch_pop_af_vcf_with_tbi // channel: [ val(meta), path(vcf.gz), path(vcf.gz.tbi) ]
    ch_svc_models // channel: [ val(meta), path(model_snv), path(model_indel) ] or [ [], [], [] ] for empty
    ch_truth_vcf_dir // channel: [ val(meta), path(dir) ] or [ [], [] ] for empty
    ch_train_bed_dir // channel: [ val(meta), path(dir) ] or [ [], [] ] for empty
    ch_test_bed_dir // channel: [ val(meta), path(dir) ] or [ [], [] ] for empty
    ch_intervals // channel: [ val(meta), path(bed) ] or [ [], [] ] for empty
    use_gatk_haplotype_caller // value: true or false
    aligner // value: aligner name (e.g. 'pb_giraffe', 'bwa')
    enable_small_variant_caller // value: true or false
    enable_small_variant_caller_retrain // value: true or false
    enable_partitioning // val: boolean
    partition_count // val: integer

    main:

    if (use_gatk_haplotype_caller) {
        HAPLOTYPECALLER_GATK(
            ch_bams,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_intervals,
            enable_partitioning,
            partition_count,
        )
        hc_vcfs = HAPLOTYPECALLER_GATK.out.vcf
        hc_bams = HAPLOTYPECALLER_GATK.out.bam
    }
    else {
        HAPLOTYPECALLER_PARABRICKS(
            ch_bams,
            ch_reference,
            ch_reference_dict,
            ch_splitintervals_reference,
            ch_splitintervals_reference_dict,
            ch_intervals,
            enable_partitioning,
            partition_count,
        )
        hc_vcfs = HAPLOTYPECALLER_PARABRICKS.out.vcf
        hc_bams = HAPLOTYPECALLER_PARABRICKS.out.bam
    }
    ch_svc_input = hc_bams.join(hc_vcfs)
    if (enable_small_variant_caller) {
        XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS(
            ch_svc_input,
            ch_reference,
            ch_svc_models,
            ch_pop_af_vcf_with_tbi,
            ch_par_beds,
            ch_intervals,
            aligner,
        )
        ch_pretrained_filtered_vcfs = XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS.out.vcf

    }
    if (enable_small_variant_caller_retrain) {
        // Match sample with sample-specific training resources (truth VCF, training BED, test BED)
        ch_truth_vcf = ch_bams
            .combine(ch_truth_vcf_dir)
            .map { bam_meta, _bam_path, _bam_index, _truth_meta, truth_dir ->
                def sample_id = bam_meta.id.split('_')[0]
                def fileList = truth_dir ? file(truth_dir).listFiles() : null
                def vcf_matching_bam = fileList != null
                    ? fileList.find { it -> !it.name.endsWith('.vcf.gz.tbi') && it.name.endsWith('.vcf.gz') && it.name.contains(sample_id) }
                    : (truth_dir ? truth_dir.resolve("${sample_id}_truth.vcf.gz") : null)
                def tbi_matching_bam = fileList != null
                    ? fileList.find { it -> it.name.endsWith('.vcf.gz.tbi') && it.name.contains(sample_id) }
                    : (truth_dir ? truth_dir.resolve("${sample_id}_truth.vcf.gz.tbi") : null)
                return [bam_meta, vcf_matching_bam, tbi_matching_bam]
            }
            .filter { it -> it[1] != null && it[2] != null }
        // filter out samples without matching truth VCF + tbi
        ch_train_bed = ch_bams
            .combine(ch_train_bed_dir)
            .map { bam_meta, _bam_path, _bam_index, _bed_meta, bed_dir ->
                def sample_id = bam_meta.id.split('_')[0]
                def fileList = bed_dir ? file(bed_dir).listFiles() : null
                def bed_matching_bam = fileList != null
                    ? fileList.find { it -> it.name.endsWith('.bed') && it.name.contains(sample_id) }
                    : (bed_dir ? bed_dir.resolve("${sample_id}_train.bed") : null)
                return [bam_meta, bed_matching_bam]
            }
            .filter { it -> it[1] != null }
        // filter out samples without matching training BED
        ch_test_bed = ch_bams
            .combine(ch_test_bed_dir)
            .map { bam_meta, _bam_path, _bam_index, _bed_meta, bed_dir ->
                def sample_id = bam_meta.id.split('_')[0]
                def fileList = bed_dir ? file(bed_dir).listFiles() : null
                def bed_matching_bam = fileList != null
                    ? fileList.find { it -> it.name.endsWith('.bed') && it.name.contains(sample_id) }
                    : (bed_dir ? bed_dir.resolve("${sample_id}_test.bed") : null)
                return [bam_meta, bed_matching_bam]
            }
            .filter { it -> it[1] != null }
        // filter out samples without matching test BED

        // Retrain only samples that have all benchmark resources. The strict
        // joins inside RETRAIN_GERMLINE_SMALL_VARIANT_MODEL abort the pipeline
        // when a sample is missing one; inner-joining here drops such samples
        // (e.g. clinical samples with no matching truth/BED) instead.
        ch_retrain_input = ch_svc_input
            .join(ch_train_bed)
            .join(ch_test_bed)
            .join(ch_truth_vcf)
            .map { meta, bams, bam_indexes, vcf, vcf_index, _train_bed, _test_bed, _truth_vcf, _truth_tbi ->
                [meta, bams, bam_indexes, vcf, vcf_index]
            }

        // Warn about dropped samples so a silent skip does not hide a misnamed
        // sample or a missing resource file.
        ch_svc_input
            .join(ch_retrain_input.map { meta, _bams, _bam_indexes, _vcf, _vcf_index -> [meta, true] }, remainder: true)
            .filter { it -> it[-1] == null }
            .subscribe { it ->
                def meta = it[0]
                log.warn("Skipping germline SVC retraining for sample '${meta.id}': no matching benchmark resources (truth VCF, training BED, test BED) for sample prefix '${meta.id.split('_')[0]}'.")
            }

        RETRAIN_GERMLINE_SMALL_VARIANT_MODEL(
            ch_retrain_input,
            ch_reference,
            ch_par_beds,
            ch_pop_af_vcf_with_tbi,
            ch_train_bed,
            ch_test_bed,
            ch_truth_vcf,
            aligner,
        )
        ch_retrained_vcfs = RETRAIN_GERMLINE_SMALL_VARIANT_MODEL.out.vcf
        ch_retrained_models = RETRAIN_GERMLINE_SMALL_VARIANT_MODEL.out.models
    }

    emit:
    raw_vcf              = hc_vcfs // channel: [ val(meta), path(vcf), path(tbi) ]
    vcf_pretrained_model = enable_small_variant_caller ? ch_pretrained_filtered_vcfs : channel.empty() // channel: [ val(meta), path(vcf), path(tbi) ]
    vcf_retrained_model  = enable_small_variant_caller_retrain ? ch_retrained_vcfs : channel.empty() // channel: [ val(meta), path(vcf), path(tbi) ]
    metrics_pretrained_model = enable_small_variant_caller ? XOOS_SMALL_VARIANT_CALLER_FILTER_VARIANTS_GERMLINE_WGS.out.metrics : channel.empty() // channel: [ val(meta), path(metrics_tsv) ]
    metrics_retrained_model  = enable_small_variant_caller_retrain ? RETRAIN_GERMLINE_SMALL_VARIANT_MODEL.out.metrics : channel.empty() // channel: [ val(meta), path(metrics_tsv) ]
    retrained_models     = enable_small_variant_caller_retrain ? ch_retrained_models : channel.empty() // channel: [ val(meta), path(snv), path(indel) ]
}
