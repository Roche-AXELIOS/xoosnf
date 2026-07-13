process ADD_RUN_INFO_TO_MULTIQC_HEADER {
    conda false
    container 'ghcr.io/roche-axelios/xoosnf/add-run-info-to-multiqc-header:1.1.0'

    label 'process_single'

    tag "${meta.id}"

    input:
    tuple val(meta), path(multiqc_config)

    output:
    tuple val(meta), path("multiqc_config_with_run_info.yml"), emit: multiqc_config_with_run_info
    tuple val("${task.process}"), val('add-run-info-to-multiqc-header'), eval('add-run-info-to-multiqc-header --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    add-run-info-to-multiqc-header \\
        --input-config ${multiqc_config} \\
        --output-config multiqc_config_with_run_info.yml \\
        --run-name '${meta.id}' \\
        --run-type '${meta.run_type}'
    """

    stub:
    """
    touch multiqc_config_with_run_info.yml
    """
}
