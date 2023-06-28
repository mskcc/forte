include { HTSEQ_COUNT    } from '../../modules/local/htseq/count/main'
include { KALLISTO_QUANT } from '../../modules/local/kallisto/quant/main'


workflow QUANTIFICATION {
    take:
    bam
    bai
    gtf
    reads
    kallisto_idx

    main:

    ch_versions = Channel.empty()

    HTSEQ_COUNT(
        bam.join(bai,by:[0]),
        gtf
    )
    ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions)

    KALLISTO_QUANT(
        reads,
        kallisto_idx,
        [],
        []
    )
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)

    emit:
    htseq_counts    = HTSEQ_COUNT.out.counts
    ch_versions     = ch_versions
}
