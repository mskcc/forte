include { STAR_ALIGN       } from '../../modules/nf-core/star/align/main'
include { SAMTOOLS_MERGE   } from '../../modules/nf-core/samtools/merge/main'
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
        false,
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    star_align_bam = STAR_ALIGN.out.bam
        .map{ meta, bam ->
	    def meta_clone = meta.clone()
	    meta_clone.id = meta.sample
	    [meta_clone, bam]
	}.branch { meta, bam ->
            needs_merge: meta.fq_num > 1
	    skips_merge: meta.fq_num == 1
	}

    SAMTOOLS_MERGE(
        star_align_bam.needs_merge
            .groupTuple(by: [0]),
        [],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    merged_bam = star_align_bam.skips_merge
	.mix(SAMTOOLS_MERGE.out.bam)
    
    SAMTOOLS_INDEX(merged_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    UMITOOLS_DEDUP(
        merged_bam
            .filter{ meta, bam -> meta.has_umi }
            .join(SAMTOOLS_INDEX.out.bai, by:[0]),
        true
    )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

    SAMTOOLS_INDEX_DEDUP(UMITOOLS_DEDUP.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions.first())

    dedup_bam = UMITOOLS_DEDUP.out.bam 
        .mix(
            merged_bam
                .filter{ meta, bam -> ! meta.has_umi }
        )


    emit:
    bam             = dedup_bam 
    bam_withdup     = merged_bam
    bai             = SAMTOOLS_INDEX_DEDUP.out.bai
    bai_withdup     = SAMTOOLS_INDEX.out.bai
    ch_versions     = ch_versions
}
