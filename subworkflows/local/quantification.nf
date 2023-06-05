include { HTSEQ_COUNT } from '../../modules/local/htseq/count/main'

workflow QUANTIFICATION {
    take:
    bam
    bai
    gtf

    main:

    ch_versions = Channel.empty()

    HTSEQ_COUNT(bam.join(bai,by:[0]),gtf)
    ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions)

    emit:
    htseq_counts = HTSEQ_COUNT.out.counts
    ch_versions = ch_versions
}
