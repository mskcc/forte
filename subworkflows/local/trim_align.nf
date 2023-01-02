include { FASTP            } from '../../modules/nf-core/fastp/main'
include { STAR_ALIGN       } from '../../modules/nf-core/star/align/main'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/umitools/extract/main'
include { UMITOOLS_DEDUP   } from '../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX   } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_BAM2FQ  } from '../../modules/nf-core/samtools/bam2fq/main'

workflow TRIM_ALIGN {

    take:
    reads
    star_index
    gtf
    run_umitools

    main:

    ch_versions = Channel.empty()

    adapter_fasta = file("${workflow.projectDir}/assets/adapters.fasta")
    if ( ! adapter_fasta.exists() ){ adapter_fasta = [] }

    if (run_umitools){
        UMITOOLS_EXTRACT(reads)
        reads4fastp = UMITOOLS_EXTRACT.out.reads
    } else {
        reads4fastp = reads
    }

    FASTP(reads4fastp, adapter_fasta, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    if (run_umitools || params.run_alignment) {
    STAR_ALIGN(
        FASTP.out.reads,
        star_index,
        gtf,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    processed_bais = SAMTOOLS_INDEX.out.bai
    } else {
        processed_bams = Channel.empty()
        processed_bais = Channel.empty()
    }

    if (run_umitools){
        UMITOOLS_DEDUP(
            STAR_ALIGN.out.bam
                .join(SAMTOOLS_INDEX.out.bai, by:[0]),
                false
        )
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

        SAMTOOLS_BAM2FQ(UMITOOLS_DEDUP.out.bam,true)
        ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())
        processed_reads = SAMTOOLS_BAM2FQ.out.reads
            .map{ meta, reads ->
                [meta, reads.findAll{ !(it.getName().endsWith("singleton.fq.gz") || it.getName().endsWith("other.fq.gz")) }]
            }
        processed_bams = UMITOOLS_DEDUP.out.bam
    } else {
        processed_reads = FASTP.out.reads
    }

    emit:
    bam             = processed_bams
    bai             = processed_bais
    reads           = processed_reads
    fastp_json      = FASTP.out.json
    ch_versions     = ch_versions

}
