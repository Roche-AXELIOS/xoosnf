include { KMC_KMER_COUNT } from '../../../modules/local/kmc/main'
include { VG_HAPLOTYPES  } from '../../../modules/local/vg/haplotypes/main'
include { VG_AUTOINDEX   } from '../../../modules/local/vg/autoindex/main'

workflow CREATE_PERSONALIZED_PANGENOME {
    take:
    ch_fastq // channel: [ val(meta), path(reads) ]
    ch_giraffe_index // channel: [ [id:'giraffe_idx'], path(idx_dir) ]

    main:
    // If ch_fastq contains the field `personalized_on` in the meta
    // we use the path specified in that field as the input to KMC, otherwise we use the fastq path.
    // `personalized_on` is already a list of staged file paths (assigned from `fastqs` in xoos.nf),
    // so we pass it through directly rather than wrapping it in file().
    ch_fastq_for_personalization = ch_fastq.map { meta, reads ->
        if (meta.containsKey('personalized_on') && meta.personalized_on) {
            return [meta, meta.personalized_on]
        }
        else {
            return [meta, reads]
        }
    }

    kmc_out = KMC_KMER_COUNT(ch_fastq_for_personalization, params.pangenome_personalization_downsample_read_count)

    haplotypes = VG_HAPLOTYPES(kmc_out.kff, ch_giraffe_index).graphs
    graphs = VG_AUTOINDEX(haplotypes).graphs

    emit:
    graphs = graphs // channel: [ val(meta), path(gbz), path(zipcodes), path(minimizer), path(dist) ]
}
