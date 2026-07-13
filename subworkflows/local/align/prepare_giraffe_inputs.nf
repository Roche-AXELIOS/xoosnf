include { CREATE_PERSONALIZED_PANGENOME } from './create_personalized_pangenome'

def resolveGiraffeIndexComponents(idxDir) {
    // Use file() with type:'dir' so Nextflow resolves the path through
    // the correct filesystem provider (local or S3/GCS).  Path.list()
    // delegates to that provider, unlike File.listFiles() which only
    // works on local filesystems.
    def baseDir = file(idxDir, type: 'dir', checkIfExists: true)
    def gbzFiles = baseDir.list()?.findAll { it -> it.endsWith('.gbz') }?.sort()

    if (!gbzFiles) {
        error("No .gbz file found in Giraffe index directory: ${idxDir}")
    }

    def prefix = gbzFiles.first().replaceFirst(/\.gbz$/, '')

    return [
        gbz : file(baseDir.resolve("${prefix}.gbz"), checkIfExists: true),
        zipcodes : file(baseDir.resolve("${prefix}.zipcodes"), checkIfExists: true),
        minimizer : file(baseDir.resolve("${prefix}.min"), checkIfExists: true),
        dist : file(baseDir.resolve("${prefix}.dist"), checkIfExists: true),
    ]
}

workflow PREPARE_GIRAFFE_INPUTS {
    take:
    ch_fastq // [ val(meta), path(reads) ]
    ch_giraffe_index // [ val(meta), path(idx_dir) ]
    personalization_enabled // val(boolean)

    main:
    if (personalization_enabled) {
        personalized_graphs = CREATE_PERSONALIZED_PANGENOME(
            ch_fastq,
            ch_giraffe_index,
        ).graphs
        align_inputs = ch_fastq
            .join(personalized_graphs)
            .map { meta, reads, graph, zipcodes, minimizers, dist -> [meta, reads, graph, zipcodes, minimizers, dist] }
    }
    else {
        static_graphs = ch_giraffe_index
            .map { meta, idx_dir ->
                def components = resolveGiraffeIndexComponents(idx_dir)
                [meta, components.gbz, components.zipcodes, components.minimizer, components.dist]
            }

        align_inputs = ch_fastq
            .combine(static_graphs)
            .map { meta, reads, _idx_meta, graph, zipcodes, minimizers, dist -> [meta, reads, graph, zipcodes, minimizers, dist] }
    }

    emit:
    align_inputs = align_inputs // [ val(meta), path(reads), path(gbz), path(zipcodes), path(min), path(dist) ]
}
