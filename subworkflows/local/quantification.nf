include { HTSEQ_COUNT           } from '../../modules/local/htseq/count/main'
include { KALLISTO_QUANT        } from '../../modules/nf-core/kallisto/quant/main'
include { SUBREAD_FEATURECOUNTS } from '../../modules/nf-core/subread/featurecounts/main'
include { COUNT_FEATURES           } from '../../modules/local/count_features/main'


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
    ch_versions   = ch_versions.mix(HTSEQ_COUNT.out.versions)


    SUBREAD_FEATURECOUNTS(
        bam,
        gtf
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    KALLISTO_QUANT(
        reads,
        kallisto_idx,
        [],
        []
    )
    ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)

    COUNT_FEATURES(
        KALLISTO_QUANT.out.abundance,
        gtf
    )


    emit:
    htseq_counts           = HTSEQ_COUNT.out.counts
    htseq_summary          = HTSEQ_COUNT.out.summary
    featurecounts_counts   = SUBREAD_FEATURECOUNTS.out.counts
    featurecounts_summary  = SUBREAD_FEATURECOUNTS.out.summary
    kallisto_log           = KALLISTO_QUANT.out.log
    kallisto_count_feature = COUNT_FEATURES.out.kallisto_count_feature
    ch_versions            = ch_versions
}
