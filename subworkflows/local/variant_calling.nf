include { PICARD_MARKDUPLICATES                   } from '../../modules/nf-core/picard/markduplicates/main'
include { GATK4_SPLITNCIGARREADS                  } from '../../modules/nf-core/gatk4/splitncigarreads/main'                                                                                                                                               
include { SAMTOOLS_INDEX as  SAMTOOLS_INDEX_CIGAR } from '../../modules/nf-core/samtools/index/main'                                                                                                                                               
include { SAMTOOLS_INDEX as  SAMTOOLS_INDEX_BQSR  } from '../../modules/nf-core/samtools/index/main'                                                                                                                                               
include { GATK4_BASERECALIBRATOR                  } from '../../modules/nf-core/gatk4/baserecalibrator/main'                                                                                                                                               
include { GATK4_APPLYBQSR                         } from '../../modules/nf-core/gatk4/applybqsr/main'                                                                                                                                                             
include { GATK4_MUTECT2                           } from '../../modules/nf-core/gatk4/mutect2/main'                                                                                                                                                                 

workflow VARIANT_CALLING {

    take:
        bam
        fasta
        fai
        dict

    main:
    ch_versions = Channel.empty()
    
    PICARD_MARKDUPLICATES(
        bam,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())


    GATK4_SPLITNCIGARREADS(
        PICARD_MARKDUPLICATES.out.bam
            .join(PICARD_MARKDUPLICATES.out.bai,by:0)
            .map{ meta, bam, bai ->
                [ meta, bam, bai, []]
            },
        fasta,
        fai,
        dict
    )
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())


    SAMTOOLS_INDEX_CIGAR(
        GATK4_SPLITNCIGARREADS.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_CIGAR.out.versions.first())

    GATK4_BASERECALIBRATOR(
        GATK4_SPLITNCIGARREADS.out.bam
            .join(SAMTOOLS_INDEX_CIGAR.out.bai, by:0)
            .map{ meta, bam, bai ->
                [meta, bam, bai, [] ]
            },
        fasta,
        fai,
        dict,
        [[id:params.genome],params.dbsnp],
        [[id:params.genome],params.dbsnpIndex]
    )
    ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())

    GATK4_APPLYBQSR(
        GATK4_SPLITNCIGARREADS.out.bam
            .join(SAMTOOLS_INDEX_CIGAR.out.bai, by:0)
            .join(GATK4_BASERECALIBRATOR.out.table, by:0)
            .map{ meta, bam, bai, table ->
                [ meta, bam, bai , table, [] ]
            },
        fasta
            .map{ meta, fasta -> 
            return fasta 
            },
        fai
            .map{ meta, fai -> 
                return fai 
            },
        dict
            .map{ meta, dict -> 
                return dict 
            }
    )
    ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

    SAMTOOLS_INDEX_BQSR(
        GATK4_APPLYBQSR.out.bam
    )

    GATK4_MUTECT2(
        GATK4_APPLYBQSR.out.bam
            .join(SAMTOOLS_INDEX_BQSR.out.bai, by:0)
            .map{ meta, bam, bai ->
                [meta, bam, bai, [] ]
            },
        fasta,
        fai,
        dict,
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    emit:
    ch_versions     = ch_versions
}