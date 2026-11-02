include { RTGTOOLS_FORMAT             } from '../../../modules/local/rtgtools/format/main'
include { RTGTOOLS_VCFEVAL            } from '../../../modules/local/rtgtools/vcfeval/main'
include { XOOS_PARSE_VCF_EVAL_METRICS } from '../../../modules/local/xoos_generate_metrics_report/parse_vcf_eval_metrics/main'
include { XOOS_VCFEVAL_METRICS_TO_MULTIQC  } from '../../../modules/local/xoos_metrics_to_multiqc/vcfeval_metrics/main'
include { resolveTruthSampleName      } from '../utils_samplesheet/main'

workflow EVALUATE_SMALL_VARIANT_CALLS {
    take:
    ch_vcf_tbi // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_reference // channel: [ val(meta), path(fasta), path(fai) ]
    ch_truth_vcf_dir // channel: [ val(meta), path(dir) ] or [ [], [] ] for empty
    ch_eval_bed_dir // channel: [ val(meta), path(dir) ] or [ [], [] ] for empty
    reference_sdf // val: path string to pre-built SDF directory, or null

    main:
    // Collect versions (also published via topic channels)
    ch_versions = Channel.empty()

    // Use pre-built SDF if provided, otherwise convert reference FASTA on the fly
    if (reference_sdf) {
        sdf = channel.value([[id: 'reference_sdf'], file(reference_sdf, type: 'dir', checkIfExists: true)])
    } else {
        ch_ref_for_sdf = ch_reference.map { _meta, fasta, _fai -> [[id: fasta.baseName], fasta, [], []] }
        sdf = RTGTOOLS_FORMAT(ch_ref_for_sdf).sdf
    }

    ch_truth_vcf = ch_vcf_tbi
        .combine(ch_truth_vcf_dir)
        .map { vcf_meta, _vcf, _tbi, _truth_meta, truth_dir ->
            def truth_sample_name = resolveTruthSampleName(vcf_meta)
            def vcf_matching_bam = truth_dir
                ? file(truth_dir).listFiles().find { it -> it.name.endsWith('.vcf.gz') && it.name.contains(truth_sample_name) }
                : null
            def tbi_matching_bam = truth_dir
                ? file(truth_dir).listFiles().find { it -> it.name.endsWith('.vcf.gz.tbi') && it.name.contains(truth_sample_name) }
                : null
            return [vcf_meta, vcf_matching_bam, tbi_matching_bam]
        }
        .filter { it -> it[1] != null && it[2] != null }
    // filter out samples without matching truth VCF + tbi
    ch_eval_bed_dir = ch_vcf_tbi
        .combine(ch_eval_bed_dir)
        .map { vcf_meta, _vcf, _tbi, _bed_meta, bed_dir ->
            def truth_sample_name = resolveTruthSampleName(vcf_meta)
            def bed_matching_vcf = bed_dir
                ? file(bed_dir).listFiles().find { it -> it.name.endsWith('.bed') && it.name.contains(truth_sample_name) }
                : null
            return [vcf_meta, bed_matching_vcf]
        }
        .filter { it -> it[1] != null }
    // filter out samples without matching test BED

    // Prepare truth VCF and evaluation regions BED from truth bundle
    ch_vcfeval_input = ch_vcf_tbi
        .join(ch_truth_vcf)
        .join(ch_eval_bed_dir)
        .map { meta, vcf, tbi, truth_vcf, truth_tbi, eval_bed ->
            [meta, vcf, tbi, truth_vcf, truth_tbi, eval_bed, []]
        }

    // Run RTG vcfeval
    ch_vcfeval_results = RTGTOOLS_VCFEVAL(ch_vcfeval_input, sdf)

    // Generate reports from vcfeval output directory (a real staged path, not a string-derived parent)
    vcf_eval_report = XOOS_PARSE_VCF_EVAL_METRICS(
        ch_vcfeval_results.vcfeval_outdir.combine(channel.of('snp', 'non_snp')).map { meta, vcfeval_dir, variant_type ->
            [meta, variant_type, vcfeval_dir]
        }
    ).report

    ch_multiqc_files = XOOS_VCFEVAL_METRICS_TO_MULTIQC(
        ch_vcfeval_results.summary.join(ch_vcfeval_results.snp_roc).join(ch_vcfeval_results.non_snp_roc)
    ).multiqc_report_dir

    ch_versions = ch_versions.mix(XOOS_PARSE_VCF_EVAL_METRICS.out.versions)
    ch_versions = ch_versions.mix(XOOS_VCFEVAL_METRICS_TO_MULTIQC.out.versions)

    emit:
    report               = vcf_eval_report // channel: [ val(meta), path(report) ]
    custom_multiqc_files = ch_multiqc_files // channel: [ val(meta), path(multiqc_report_dir) ]
    vcfeval_output       = ch_vcfeval_results.vcfeval_outdir // channel: [ val(meta), path(vcfeval_dir) ]
    versions             = ch_versions
}
