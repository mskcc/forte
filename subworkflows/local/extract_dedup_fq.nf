include { SAMTOOLS_BAM2FQ } from '../../modules/nf-core/samtools/bam2fq/main'

workflow EXTRACT_DEDUP_FQ {
    take:
    bam

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_BAM2FQ(
        bam,
        true
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    dedup_reads = SAMTOOLS_BAM2FQ.out.reads
        .map{ meta, reads ->
            [meta, reads.findAll{ !(it.getName().endsWith("singleton.fq.gz") || it.getName().endsWith("other.fq.gz")) }]
        }

    emit:
    dedup_reads  = dedup_reads
    ch_versions  = ch_versions
}
