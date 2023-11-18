include { STAR_ALIGN       } from '../../modules/nf-core/star/align/main'
include { UMITOOLS_DEDUP   } from '../../modules/nf-core/umitools/dedup/main'
include {
    SAMTOOLS_INDEX;
    SAMTOOLS_INDEX as SAMTOOLS_INDEX_DEDUP
} from '../../modules/nf-core/samtools/index/main'
include { INFER_STRAND     } from './infer_strand'
include {
    AMEND_STRAND as AMEND_STRAND_BAM ;
    AMEND_STRAND as AMEND_STRAND_BAI ;
    AMEND_STRAND as AMEND_STRAND_STAR_LOG
} from './infer_strand'

workflow ALIGN_READS {

    take:
    reads
    star_index
    gtf
    refflat
    fasta


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

    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    INFER_STRAND(
        STAR_ALIGN.out.bam.filter{ meta, bam ->
            meta.strandedness == "auto"
        },
        SAMTOOLS_INDEX.out.bai.filter{ meta, bai ->
            meta.strandedness == "auto"
        },
        refflat,
        fasta
    )
    ch_versions = ch_versions.mix(INFER_STRAND.out.ch_versions)

    AMEND_STRAND_BAM(
        STAR_ALIGN.out.bam,
        INFER_STRAND.out.inferred_strand
    )
    star_bam = AMEND_STRAND_BAM.out.amended_ch

    AMEND_STRAND_BAI(
        SAMTOOLS_INDEX.out.bai,
        INFER_STRAND.out.inferred_strand
    )
    star_bai = AMEND_STRAND_BAI.out.amended_ch

    AMEND_STRAND_STAR_LOG(
        STAR_ALIGN.out.log_final
            .map{meta, log_final ->
                def meta_clone = meta.clone().findAll { !["read_group","fastq_pair_id"].contains(it.key) }
                [meta_clone,log_final]
            },
        INFER_STRAND.out.inferred_strand
    )
    star_log_final = AMEND_STRAND_STAR_LOG.out.amended_ch

    UMITOOLS_DEDUP(
        STAR_ALIGN.out.bam
            .filter{ meta, bam -> meta.has_umi }
            .join(SAMTOOLS_INDEX.out.bai, by:[0]),
        true
    )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

    SAMTOOLS_INDEX_DEDUP(UMITOOLS_DEDUP.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_DEDUP.out.versions.first())

    final_bam = UMITOOLS_DEDUP.out.bam
        .mix(
            STAR_ALIGN.out.bam
                .filter{ meta, bam -> ! meta.has_umi }
        )
    final_bai = SAMTOOLS_INDEX_DEDUP.out.bai
        .mix(
            SAMTOOLS_INDEX.out.bai
                .filter{ meta, bai -> ! meta.has_umi }
        )

    emit:
    bam             = final_bam
    bam_dedup       = UMITOOLS_DEDUP.out.bam
    bam_withdup     = star_bam
    bai             = final_bai
    bai_dedup       = SAMTOOLS_INDEX_DEDUP.out.bai
    bai_withdup     = star_bai
    star_log_final  = star_log_final
    umitools_dedup_log = UMITOOLS_DEDUP.out.log
    inferred_strand = INFER_STRAND.out.inferred_strand
    ch_versions     = ch_versions
}
