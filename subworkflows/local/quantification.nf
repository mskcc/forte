include { HTSEQ_COUNT } from '../../modules/local/htseq/count/main'

workflow QUANTIFICATION {
    take:
    bam
    gtf

    main:
    
    ch_versions = Channel.empty()

    HTSEQ_COUNT(bam,gtf)
    ch_versions = ch_versions.mix(HTSEQ_COUNT.out.versions)

    emit:
    ch_versions = ch_versions
}