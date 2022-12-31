include { FASTP            } from '../../modules/nf-core/fastp/main'
include { STAR_ALIGN       } from '../../modules/nf-core/star/align/main'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/umitools/extract/main'
include { UMITOOLS_DEDUP   } from '../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX   } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_BAM2FQ  } from '../../modules/nf-core/samtools/bam2fq/main'

workflow TRIM_ALIGN_WITH_UMI {

    take:
    reads
    star_index
    gtf

    main:

    ch_versions = Channel.empty()

    adapter_fasta = file("${workflow.projectDir}/assets/adapters.fasta")
    if ( ! adapter_fasta.exists() ){ adapter_fasta = [] }
    
    UMITOOLS_EXTRACT(reads)

    FASTP(UMITOOLS_EXTRACT.out.reads, adapter_fasta, false, false)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    
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

    UMITOOLS_DEDUP(
        STAR_ALIGN.out.bam
            .join(SAMTOOLS_INDEX.out.bai, by:[0]),
            false
    )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions.first())

    SAMTOOLS_BAM2FQ(UMITOOLS_DEDUP.out.bam,true)
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())
    
	SAMTOOLS_BAM2FQ.out.reads
	    .map{ meta, reads ->
	        [meta, reads.findAll{ !(it.getName().endsWith("singleton.fq.gz") || it.getName().endsWith("other.fq.gz")) }]
		}
    ).set{processed_reads}

    processed_bams = STAR_ALIGN.out.bam
	                      .filter{ meta, bam -> ! meta.has_umi }
						  .mix(UMITOOLS_DEDUP.out.bam)

    emit:
    bam             = UMITOOLS_DEDUP.out.bam
    bai             = SAMTOOLS_INDEX.out.bai
    processed_reads = processed_reads
	fastp_json      = FASTP.out.json
    ch_versions     = ch_versions

}
