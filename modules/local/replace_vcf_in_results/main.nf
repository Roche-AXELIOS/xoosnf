process REPLACE_VCF_IN_RESULTS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(results_dir, stageAs: 'input_results/*'), path(reheadered_vcf, stageAs: 'reheadered.vcf.gz'), path(reheadered_index, stageAs: 'reheadered.vcf.gz.tbi'), val(vcf_pattern), val(result_dir_name)

    output:
    tuple val(meta), path("${result_dir_name}/*"), emit: results_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    mkdir -p '${result_dir_name}'
    cp -a input_results/. '${result_dir_name}'/
    cp reheadered.vcf.gz '${result_dir_name}'/${vcf_pattern}
    cp reheadered.vcf.gz.tbi '${result_dir_name}'/${vcf_pattern}.tbi
    """

    stub:
    """
    mkdir -p '${result_dir_name}'
    cp -a input_results/. '${result_dir_name}'/
    cp reheadered.vcf.gz '${result_dir_name}'/${vcf_pattern}
    cp reheadered.vcf.gz.tbi '${result_dir_name}'/${vcf_pattern}.tbi
    """
}
