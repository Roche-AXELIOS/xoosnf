include { CREATE_PERSONALIZED_PANGENOME } from './create_personalized_pangenome'
include { resolveGiraffeIndexComponents } from '../utils_xoos_workflow/main'

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
