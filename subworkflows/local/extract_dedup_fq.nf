include { GATK4_SAMTOFASTQ } from '../../modules/nf-core/gatk4/samtofastq/main'   

workflow EXTRACT_DEDUP_FQ {
    take:
    bam

    main:
    ch_versions = Channel.empty()

    GATK4_SAMTOFASTQ(
        bam
    )
    ch_versions = ch_versions.mix(GATK4_SAMTOFASTQ.out.versions.first())

    dedup_reads = GATK4_SAMTOFASTQ.out.fastq

    emit:
    dedup_reads  = dedup_reads
    ch_versions  = ch_versions
}
