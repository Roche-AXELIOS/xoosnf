process UNTAR {
    tag "$meta.id"

    input:
    tuple val(meta), path(archive)

    output:
    tuple val(meta), path("${prefix}"), emit: untar

    script:
    prefix = archive.baseName.replaceAll(/\.tar$/, '')
    """
    mkdir -p ${prefix}
    tar xzf ${archive} -C ${prefix} --strip-components 1
    """
}
