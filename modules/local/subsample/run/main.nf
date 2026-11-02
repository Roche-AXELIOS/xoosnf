process SUBSAMPLE_RUN {
    tag "${meta.id}"
    label 'process_single'

    conda false
    container 'ghcr.io/roche-axelios/xoosnf/subsample-run:1.2.0'

    input:
    tuple val(meta), path(input_dir), path(samplesheet)

    output:
    tuple val(meta), path("subsampled/${input_dir}"), path(samplesheet), emit: subsampled
    tuple val("${task.process}"), val('subsample-run'), eval('subsample-run --version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def real_input_dir = input_dir.toRealPath()
    """
    mkdir -p 'subsampled/${input_dir}'

    subsample-run \\
        --input-dir '${real_input_dir}' \\
        --output-dir 'subsampled/${input_dir}' \\
        ${args}
    """

    stub:
    """
    mkdir -p 'subsampled/${input_dir}'
    touch 'subsampled/${input_dir}/stub.txt'
    """
}
