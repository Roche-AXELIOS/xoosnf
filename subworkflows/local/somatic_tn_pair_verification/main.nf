include { SOMALIER_EXTRACT as SOMALIER_EXTRACT_TUMOR; SOMALIER_EXTRACT as SOMALIER_EXTRACT_NORMAL } from '../../../modules/nf-core/somalier/extract/main.nf'
include { SOMALIER_RELATE                              } from '../../../modules/nf-core/somalier/relate/main.nf'
include { XOOS_SOMALIER_METRICS_TO_MULTIQC             } from '../../../modules/local/xoos_metrics_to_multiqc/somalier_metrics/main.nf'
include { XOOS_PARSE_SOMALIER_METRICS as XOOS_PARSE_SOMALIER_PAIR_METRICS   } from '../../../modules/local/xoos_generate_metrics_report/parse_somalier_metrics/main.nf'
include { XOOS_PARSE_SOMALIER_METRICS as XOOS_PARSE_SOMALIER_SAMPLE_METRICS } from '../../../modules/local/xoos_generate_metrics_report/parse_somalier_metrics/main.nf'

workflow SOMATIC_TN_PAIR_VERIFICATION {

    take:
    ch_tn_bam_pairs             // channel: [ val(meta), [tumor_bam, normal_bam], [tumor_bai, normal_bai] ]
    ch_reference                // channel: [ val(meta_ref), path(fasta), path(fai) ]
    ch_sites_vcf                // channel: [ val(meta), path(sites_vcf) ]
    somatic_tn_pair_mismatch_action     // val: 'warn', 'abort', or 'skip'
    concordance_cutoff          // val: concordance threshold (double)

    main:

    ch_versions = Channel.empty()

    // Split reference into fasta and fai channels (matching SOMALIER_EXTRACT input signature)
    ch_fasta = ch_reference.map { meta, fasta, _fai -> [meta, fasta] }
    ch_fai   = ch_reference.map { meta, _fasta, fai -> [meta, fai] }

    // Split T/N pairs into individual BAMs for per-sample extraction
    ch_tumor_bams = ch_tn_bam_pairs.map { meta, bams, bais ->
        [[id: meta.tumor_id, group_id: meta.id, run_name: meta.tumor_run_name, run_type: meta.run_type], bams[0], bais[0]]
    }
    ch_normal_bams = ch_tn_bam_pairs.map { meta, bams, bais ->
        [[id: meta.normal_id, group_id: meta.id, run_name: meta.normal_run_name, run_type: meta.run_type], bams[1], bais[1]]
    }

    // Run somalier extract on tumor and normal BAMs
    SOMALIER_EXTRACT_TUMOR(ch_tumor_bams, ch_fasta, ch_fai, ch_sites_vcf)
    SOMALIER_EXTRACT_NORMAL(ch_normal_bams, ch_fasta, ch_fai, ch_sites_vcf)

    // Group extracted .somalier files by T/N group for relate
    ch_extracts = SOMALIER_EXTRACT_TUMOR.out.extract
        .map { meta, somalier_file -> [[id: meta.group_id], [id: meta.id, run_name: meta.run_name, run_type: meta.run_type], somalier_file] }
        .combine(
            SOMALIER_EXTRACT_NORMAL.out.extract
                .map { meta, somalier_file -> [[id: meta.group_id], [id: meta.id, run_name: meta.run_name, run_type: meta.run_type], somalier_file] },
            by: 0
        )
        .map { group_key, tumor_meta, tumor_somalier, normal_meta, normal_somalier ->
            // SOMALIER_RELATE expects: [ val(meta), path(extract_files), path(ped) ]
            def group_meta = [
                id: group_key.id,
                run_name: tumor_meta.run_name,
                tumor_id: tumor_meta.id,
                normal_id: normal_meta.id,
                tumor_run_name: tumor_meta.run_name,
                normal_run_name: normal_meta.run_name,
                run_type: tumor_meta.run_type]
            [group_meta, [tumor_somalier, normal_somalier], []]
        }

    // Run somalier relate per T/N group (no PED file, no sample groups)
    SOMALIER_RELATE(ch_extracts, [])

    // Collect versions from somalier processes
    ch_versions = ch_versions.mix(SOMALIER_EXTRACT_TUMOR.out.versions_somalier.map { _process, tool, version -> "somalier_extract:\n  ${tool}: ${version}" })
    ch_versions = ch_versions.mix(SOMALIER_RELATE.out.versions_somalier.map { _process, tool, version -> "somalier_relate:\n  ${tool}: ${version}" })

    // Parse concordance from pairs.tsv and gate downstream processing.
    // Join concordance results back to the original T/N pairs so downstream
    // channels are blocked until the somalier check completes.
    ch_concordance = SOMALIER_RELATE.out.pairs_tsv
        .map { meta, pairs_tsv ->
            def lines = pairs_tsv.readLines()
            if (lines.size() > 1) {
                def header = lines[0].tokenize('\t')
                def data   = lines[1].tokenize('\t')
                def idx    = header.indexOf('hom_concordance')
                if (idx >= 0 && idx < data.size()) {
                    return [meta.id, data[idx] as Double]
                }
            }
            return [meta.id, null]
        }

    ch_verified_pairs = ch_tn_bam_pairs
        .map { meta, bams, bais -> [meta.id, meta, bams, bais] }
        .join(ch_concordance, by: 0)
        .map { group_id, meta, bams, bais, concordance ->
            if (concordance == null) {
                log.warn("Tumor-normal pair verification: Could not parse concordance for T/N pair '${group_id}'.")
            } else if (concordance < concordance_cutoff) {
                def msg = "Tumor-normal pair verification: T/N pair '${group_id}' concordance ${concordance} is below threshold ${concordance_cutoff}. Samples may not be from the same individual."
                if (somatic_tn_pair_mismatch_action == 'abort') {
                    error(msg)
                } else if (somatic_tn_pair_mismatch_action == 'skip') {
                    log.warn(msg + " Skipping pair.")
                    return null
                } else {
                    log.warn(msg)
                }
            } else {
                log.info("Tumor-normal pair verification: T/N pair '${group_id}' concordance = ${concordance} (PASS)")
            }
            [meta, bams, bais]
        }
        .filter { it -> it != null }

    ch_somalier_output = SOMALIER_RELATE.out.pairs_tsv
        .join(SOMALIER_RELATE.out.samples_tsv)
        .map { meta, pairs, samples -> [meta, pairs, samples] }

    ch_somalier_publish_output = ch_somalier_output
        .join(SOMALIER_RELATE.out.html)
        .map { meta, pairs, samples, html -> [meta, pairs, samples, html] }

    XOOS_SOMALIER_METRICS_TO_MULTIQC(ch_somalier_output)

    ch_multiqc_files = XOOS_SOMALIER_METRICS_TO_MULTIQC.out.multiqc_report_dir
        .map { meta, multiqc_report_dir -> [[id: meta.tumor_run_name, run_type: meta.run_type], multiqc_report_dir] }

    ch_pair_input = SOMALIER_RELATE.out.pairs_tsv
        .map { meta, pairs ->
            def pair_meta = [
                id:              meta.id,
                run_name:        meta.tumor_run_name,
                run_type:        meta.run_type,
            ]
            [pair_meta, pairs, [], []]
        }

    XOOS_PARSE_SOMALIER_PAIR_METRICS(ch_pair_input)

    // Parse sample metrics — one task per run.
    // A samples.tsv holds rows for both the tumor and the normal, which may
    // belong to different runs, so each file is fanned out to every run it
    // contributes a sample to. Each task is given only the sample ids for its
    // own run and filters the rows accordingly, keeping run attribution correct.
    ch_sample_input = SOMALIER_RELATE.out.samples_tsv
        .flatMap { meta, samples ->
            [
                [[id: meta.tumor_run_name, run_name: meta.tumor_run_name, run_type: meta.run_type], meta.tumor_id, samples],
                [[id: meta.normal_run_name, run_name: meta.normal_run_name, run_type: meta.run_type], meta.normal_id, samples],
            ]
        }
        .groupTuple()
        // A pair whose tumor and normal share a run contributes the same
        // samples.tsv twice; de-duplicate to avoid staging collisions.
        // toUnique() is used over unique() because the latter mutates in place.
        .map { meta, sample_ids, samples_files -> [meta, [], samples_files.toUnique(), sample_ids.toUnique()] }

    XOOS_PARSE_SOMALIER_SAMPLE_METRICS(ch_sample_input)

    ch_report = XOOS_PARSE_SOMALIER_PAIR_METRICS.out.report
        .mix(XOOS_PARSE_SOMALIER_SAMPLE_METRICS.out.report)

    emit:
    verified_pairs = ch_verified_pairs                // channel: [ val(meta), [tumor_bam, normal_bam], [tumor_bai, normal_bai] ]
    multiqc_files  = ch_multiqc_files                 // channel: [ val(meta), path(dir) ]
    report         = ch_report                        // channel: [ val(meta), path(json) ]
    somalier_output = ch_somalier_publish_output      // channel: [ val(meta), path(pairs.tsv), path(samples.tsv), path(html) ]
    versions       = ch_versions                      // channel: versions strings
}
