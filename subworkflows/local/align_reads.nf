include { STAR_ALIGN       } from '../../modules/nf-core/star/align/main'
include { UMITOOLS_DEDUP   } from '../../modules/nf-core/umitools/dedup/main'
include {
    SAMTOOLS_INDEX;
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP
} from '../../modules/nf-core/samtools/index/main'

workflow ALIGN_READS {

    take:
    reads
    star_index
    gtf

    main:

    ch_versions = Channel.empty()

    STAR_ALIGN(
        reads,
        star_index,
        gtf,
        false,
        [],
        []
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    merged_bam = STAR_ALIGN.out.bam

    SAMTOOLS_INDEX(merged_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    UMITOOLS_DEDUP(
        STAR_ALIGN.out.bam
            .filter{ meta, bam -> meta.has_umi }
            .join(SAMTOOLS_INDEX.out.bai, by:[0]),
        true
    )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

    SAMTOOLS_INDEX_DEDUP(UMITOOLS_DEDUP.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions.first())

    dedup_bam = UMITOOLS_DEDUP.out.bam
        .mix(
            STAR_ALIGN.out.bam
                .filter{ meta, bam -> ! meta.has_umi }
        )

    dedup_bai = SAMTOOLS_INDEX_DEDUP.out.bai
        .mix(
            SAMTOOLS_INDEX.out.bai
                .filter{meta, bai -> ! meta.has_umi}
        )

    emit:
    bam             = dedup_bam
    bam_withdup     = STAR_ALIGN.out.bam
    bai             = dedup_bai
    bai_withdup     = SAMTOOLS_INDEX.out.bai
    ch_versions     = ch_versions
}
