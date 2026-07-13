process XOOS_PARSE_VCF_EVAL_METRICS {
    container 'ghcr.io/roche-axelios/xoosnf/generate-metrics-report:1.1.0'
    conda false

    label 'process_single'
    tag "${meta.id}, ${meta.svc_model_type}, ${type}"

    input:
    tuple val(meta), val(type), path(vcf_eval_metrics)

    output:
    tuple val(meta), path("${meta.id}.vcf_eval_${meta.svc_model_type}_${type}.json"), emit: report
    tuple val("${task.process}"), val('generate-metrics-report'), eval('generate-metrics-report --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate-metrics-report parse-vcf-eval-metrics \\
        --sample-name '${meta.id}' \\
        --source 'vcf_eval/${meta.svc_model_type}/${type}' \\
        --input '${vcf_eval_metrics}/${meta.id}.${type}_roc.tsv.gz' \\
        --output '${meta.id}.vcf_eval_${meta.svc_model_type}_${type}.json'
    """

    stub:
    """
    echo '{}' > '${meta.id}.vcf_eval_${meta.svc_model_type}_${type}.json'
    """
}
