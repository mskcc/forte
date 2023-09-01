include { HTSEQ_COUNT                                 } from '../../modules/local/htseq/count/main'
include { KALLISTO_QUANT                              } from '../../modules/nf-core/kallisto/quant/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_GENE } from '../../modules/nf-core/subread/featurecounts/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_EXON } from '../../modules/nf-core/subread/featurecounts/main'
include { COUNT_FEATURES                              } from '../../modules/local/count_features/main'


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


    FEATURECOUNTS_GENE(
        bam,
        gtf
    )
    ch_versions = ch_versions.mix(FEATURECOUNTS_GENE.out.versions)

    FEATURECOUNTS_EXON(
        bam,
        gtf
    )
    ch_versions = ch_versions.mix(FEATURECOUNTS_EXON.out.versions)


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
    htseq_counts               = HTSEQ_COUNT.out.counts
    htseq_summary              = HTSEQ_COUNT.out.summary
    featurecounts_gene_counts  = FEATURECOUNTS_GENE.out.counts
    featurecounts_gene_summary = FEATURECOUNTS_GENE.out.summary
    featurecounts_exon_counts  = FEATURECOUNTS_EXON.out.counts
    featurecounts_exon_summary = FEATURECOUNTS_EXON.out.summary
    kallisto_log               = KALLISTO_QUANT.out.log
    kallisto_count_feature     = COUNT_FEATURES.out.kallisto_count_feature
    ch_versions                = ch_versions
}
