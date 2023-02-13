include { SAMTOOLS_BAM2FQ } from '../../modules/nf-core/samtools/bam2fq/main'
include { CAT_FASTQ       } from '../../modules/nf-core/cat/fastq/main'

workflow MERGE_READS {
    take:
    reads
    bam

    main:
    ch_versions = Channel.empty()

    reads_ch = reads
        .map{ meta, reads ->
            def meta_clone = meta.clone()
 	    meta_clone.id = meta.sample
 	    [meta_clone, reads]
   	}.branch { meta, reads ->
            needs_merge: ( meta.fq_num > 1 ) && ( ! ( meta.has_umi && params.dedup_umi_for_fusions ) )
	    needs_bam2fq: meta.has_umi && params.dedup_umi_for_fusions
   	    skips_merge: true 
	}

     bam_ch = bam
        .branch { meta, bam ->
            needs_bam2fq: meta.has_umi && params.dedup_umi_for_fusions
	    skips_bam2fq: true
	}


    CAT_FASTQ(reads_ch.needs_merge)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    SAMTOOLS_BAM2FQ(
	bam_ch.needs_bam2fq,
        true
    )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    merged_reads = reads_ch.skips_merge
        .mix(CAT_FASTQ.out.reads)
        .mix(
            SAMTOOLS_BAM2FQ.out.reads
                .map{ meta, reads ->
                    [meta, reads.findAll{ !(it.getName().endsWith("singleton.fq.gz") || it.getName().endsWith("other.fq.gz")) }]
                }
        )
    
    emit:
    merged_reads = merged_reads
    ch_versions  = ch_versions
}
