include { GATK4_SPLITNCIGARREADS } from '../../modules/nf-core/gatk4/splitncigarreads/main'                                                                                                                                               
include { GATK4_BASERECALIBRATOR } from '../../modules/nf-core/gatk4/baserecalibrator/main'                                                                                                                                               
include { GATK4_APPLYBQSR        } from '../../modules/nf-core/gatk4/applybqsr/main'                                                                                                                                                             
include { GATK4_MUTECT2          } from '../../modules/nf-core/gatk4/mutect2/main'                                                                                                                                                                 

workflow VARIANT_CALLING {

    take:
        bam
        bai
        fasta
        fai
        dict

    main:
    ch_versions = Channel.empty()
    
    GATK4_SPLITNCIGARREADS(
        bam
            .combine(bai,by:[0]),
        fasta,
        fai,
        dict
    )

    GATK4_BASERECALIBRATOR(
        GATK4_SPLITNCIGARREADS.out.bam,
        fasta,
        fai,
        dict
    )

     GATK4_APPLYBQSR(
        GATK4_SPLITNCIGARREADS.out.bam
            .combine(GATK4_BASERECALIBRATOR.out.table, by:[0])
            .map{ meta, bam, table ->
                [ meta, input, '' , bqsr_table]
            },
        fasta,
        fai,
        dict
    )

}