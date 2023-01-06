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
    run_bam2fq

    main:

    ch_versions = Channel.empty()

    adapter_fasta = file("${workflow.projectDir}/assets/adapters.fasta")
    if ( ! adapter_fasta.exists() ){ adapter_fasta = [] }

    if (run_umitools){
        UMITOOLS_EXTRACT(reads)
        ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())
        extracted_reads = UMITOOLS_EXTRACT.out.reads
    } else {
        extracted_reads = Channel.empty()
    }

    // once FASTP can run without producing reads then fix this logic.
    FASTP(
        run_umitools ? extracted_reads : reads,
        adapter_fasta,
        false,
        false
    )
    if (params.skip_trimming){
        trimmed_reads = Channel.empty()
    } else {
        trimmed_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }

    reads_for_alignment = params.skip_trimming ? (run_umitools ? extracted_reads : reads) : trimmed_reads

    STAR_ALIGN(
        reads_for_alignment,
        star_index,
        gtf,
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    SAMTOOLS_INDEX(STAR_ALIGN.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    if (run_umitools){
        UMITOOLS_DEDUP(
            STAR_ALIGN.out.bam
                .join(SAMTOOLS_INDEX.out.bai, by:[0]),
            false
        )
        ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

        dedup_bam = UMITOOLS_DEDUP.out.bam

        if (run_bam2fq) {
            SAMTOOLS_BAM2FQ(UMITOOLS_DEDUP.out.bam,true)
            ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())
            dedup_reads = SAMTOOLS_BAM2FQ.out.reads
                .map{ meta, reads ->
                    [meta, reads.findAll{ !(it.getName().endsWith("singleton.fq.gz") || it.getName().endsWith("other.fq.gz")) }]
                }
        } else {
            dedup_reads = Channel.empty()
        }

    } else {
        dedup_bam = Channel.empty()
        dedup_reads = Channel.empty()
    }

    emit:
    reads           = run_bam2fq ? dedup_reads : reads_for_alignment
    bam             = run_umitools ? dedup_bam : STAR_ALIGN.out.bam
    bai             = SAMTOOLS_INDEX.out.bai
    fastp_json      = FASTP.out.json
    ch_versions     = ch_versions
}
