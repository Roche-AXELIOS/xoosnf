#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Roche-AXELIOS/xoosnf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Roche-AXELIOS/xoosnf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { XOOS                    } from './workflows/xoosnf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_xoos_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_xoos_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version.toBoolean(),
        params.validate_params.toBoolean(),
        params.monochrome_logs.toBoolean(),
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full.toBoolean(),
        params.show_hidden.toBoolean()
    )

    //
    // WORKFLOW: Run main workflow
    //
    XOOS (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email.toBoolean(),
        params.outdir,
        params.monochrome_logs.toBoolean(),
        XOOS.out.multiqc_report
    )

    publish:
    bam = XOOS.out.bam
    pre_subsample_bam = XOOS.out.pre_subsample_bam
    vcf = XOOS.out.vcf
    raw_vc_vcf = XOOS.out.raw_vc_vcf
    svc_pretrained_vcf = XOOS.out.svc_pretrained_vcf
    svc_retrained_vcf = XOOS.out.svc_retrained_vcf
    vcfeval_output = XOOS.out.vcfeval_output
    demuxed_fastqs = XOOS.out.demuxed_fastqs
    demux_metrics_files = XOOS.out.demux_metrics_files
    combined_report = XOOS.out.combined_report
    combined_report_summary = XOOS.out.combined_report_summary
    multiqc_report = XOOS.out.multiqc_report
    cnc_output = XOOS.out.cnc_output
    str_output = XOOS.out.str_output
    tfe_output = XOOS.out.tfe_output
    contamination_output = XOOS.out.contamination_output
    alignment_metrics_files = XOOS.out.alignment_metrics_files
    gc_bias_metrics_files = XOOS.out.gc_bias_metrics_files
    hs_metrics_files = XOOS.out.hs_metrics_files
    read_collapser_metrics_files = XOOS.out.read_collapser_metrics_files
}

// Each output declares an index file so that tuple-based outputs are included
// in provenance reports. JSON is used for all outputs for Nextflow 25.10.x
// compatibility (CSV with header=true requires all record elements to be maps).
output {
    bam {
        path { meta, _bam_path, _bai_path -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/bam.json'
        }
    }
    pre_subsample_bam {
        path { meta, _bam_path, _bai_path -> "${meta.run_name}/${meta.id}/pre_subsample" }
        index {
            path 'index/pre_subsample_bam.json'
        }
    }
    vcf {
        path { meta, _vcf_path, _tbi_path -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/vcf.json'
        }
    }
    raw_vc_vcf {
        path { meta, _vcf_path, _tbi_path -> "${meta.run_name}/${meta.id}/${meta.raw_variant_caller}" }
        index {
            path 'index/raw_vc_vcf.json'
        }
    }
    svc_pretrained_vcf {
        path { meta, _vcf_path, _tbi_path -> "${meta.run_name}/${meta.id}/small_variant_caller/pretrained" }
        index {
            path 'index/svc_pretrained_vcf.json'
        }
    }
    svc_retrained_vcf {
        path { meta, _vcf_path, _tbi_path -> "${meta.run_name}/${meta.id}/small_variant_caller/retrained" }
        index {
            path 'index/svc_retrained_vcf.json'
        }
    }
    vcfeval_output {
        // The emitted directory is named 'vcfeval', so publishing to the model
        // type parent yields '<model_type>/vcfeval/'. svc_model_type is e.g.
        // 'small_variant_pretrained'/'small_variant_retrained' (germline) or
        // 'pretrained' (somatic); normalise to pretrained/retrained so the raw
        // vcfeval output sits beside the corresponding SVC VCF.
        path { meta, _vcfeval_dir ->
            def model_type = (meta.svc_model_type ?: 'pretrained').replaceFirst(/^small_variant_/, '')
            "${meta.run_name}/${meta.id}/small_variant_caller/${model_type}"
        }
        index {
            path 'index/vcfeval_output.json'
        }
    }
    demuxed_fastqs {
        path { meta, _fastq_dir -> "${meta.run_name}/${meta.id}/fastqs" }
        index {
            path 'index/demuxed_fastqs.json'
        }
    }
    demux_metrics_files {
        // Run-level output (meta.id == run_name). The emitted value is the
        // 'demux/metrics' directory, whose relative structure is preserved on
        // publish, so publishing to the run root yields '<run_name>/demux/
        // metrics/' containing run_stats.tsv, sample_stats.tsv, and the
        // read-length distribution tables.
        path { meta, _metrics_dir -> "${meta.run_name}" }
        index {
            path 'index/demux_metrics_files.json'
        }
    }
    cnc_output {
        path { meta, _cnc_out_dir -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/cnc_output.json'
        }
    }
    str_output {
        path { meta, _str_out_dir -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/str_output.json'
        }
    }
    tfe_output {
        path { meta, _tfe_out_dir -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/tfe_output.json'
        }
    }
    contamination_output {
        path { meta, _contamination_out_dir -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/contamination_output.json'
        }
    }
    alignment_metrics_files {
        // One directory per metrics stage (report_prefix) and region so that the
        // identically named coverage/accuracy/read_metrics dirs do not collide.
        path { meta, _read_metrics, _coverage, _accuracy -> "${meta.run_name}/${meta.id}/alignment_metrics/${meta.report_prefix}_${meta.metrics_region}" }
        index {
            path 'index/alignment_metrics_files.json'
        }
    }
    gc_bias_metrics_files {
        path { meta, _gc_bias_dir -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/gc_bias_metrics_files.json'
        }
    }
    hs_metrics_files {
        path { meta, _hs_metrics -> "${meta.run_name}/${meta.id}/hs_metrics" }
        index {
            path 'index/hs_metrics_files.json'
        }
    }
    read_collapser_metrics_files {
        // The emitted directory is named 'read_collapser', so publishing to the
        // sample parent yields '<run_name>/<id>/read_collapser/' containing
        // summary_stats.tsv, cluster_size_distribution_summary.tsv, and
        // cluster_size_distributions.tsv.
        path { meta, _read_collapser_metrics -> "${meta.run_name}/${meta.id}" }
        index {
            path 'index/read_collapser_metrics_files.json'
        }
    }
    combined_report {
        // meta.id here is the run name
        path { meta, _report_path -> "${meta.id}" }
        index {
            path 'index/combined_report.json'
        }
    }
    combined_report_summary {
        // meta.id here is the run name
        path { meta, _report_path -> "${meta.id}" }
        index {
            path 'index/combined_report_summary.json'
        }
    }
    multiqc_report {
        // meta.id here is the run name
        path { meta, _report_path -> "${meta.id}" }
        index {
            path 'index/multiqc_report.json'
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
