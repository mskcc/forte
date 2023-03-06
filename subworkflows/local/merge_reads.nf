include { SAMTOOLS_BAM2FQ } from '../../modules/nf-core/samtools/bam2fq/main'

workflow MERGE_READS {
    take:
    reads
    bam

    main:
    ch_versions = Channel.empty()

    reads_ch = reads
        .filter{ meta, reads -> ! ( meta.has_umi && params.dedup_umi_for_fusions) }
        
    bam_ch = bam
        .filter{ meta, bam -> meta.has_umi && params.dedup_umi_for_fusions }

    SAMTOOLS_BAM2FQ(
        bam_ch, 
        true
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    merged_reads = reads_ch
        .mix(
            SAMTOOLS_BAM2FQ.out.reads
                .map{ meta, reads ->
                    [meta, reads.findAll{ !(it.getName().endsWith("singleton.fq.gz") || it.getName().endsWith("other.fq.gz")) }]
                }
        )

    emit:
    dedup_reads  = merged_reads
    ch_versions  = ch_versions
}
