process XOOS_STR_CALLER {
    conda false
    container 'ghcr.io/roche-axelios/xoos/str_caller:1.1.0'

    label 'process_single'

    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta2), path(reference), path(reference_index)
    tuple val(meta3), path(repeat_catalog)

    output:
    tuple val(meta), path("str_caller/${prefix}_str.vcf.gz"), path("str_caller/${prefix}_str.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("str_caller/*"), emit: results_dir
    tuple val("${task.process}"), val('str_caller'), eval("str_caller --version"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    def args = task.ext.args ?: ''
    def catalog_path = repeat_catalog ? repeat_catalog : '/apps/hg38.full_repeat_catalog.json'
    """
    str_caller \\
        --reference '${reference}' \\
        --str-catalog '${catalog_path}' \\
        --bam-input '${bam}' \\
        --output-dir str_caller \\
        --output-vcf '${prefix}_str.vcf.gz' \\
        --output-json '${prefix}_str.json' \\
        --sample-name '${meta.id}' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    mkdir -p str_caller
    gzip -n < /dev/null > 'str_caller/${prefix}_str.vcf.gz'
    touch 'str_caller/${prefix}_str.vcf.gz.tbi'
    echo '{"sample": "${prefix}"}' > 'str_caller/${prefix}_str.json'
    """
}
