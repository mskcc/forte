include { GATK4_SPLITNCIGARREADS } from '../modules/nf-core/gatk4/splitncigarreads/main'                                                                                                                                               
include { GATK4_BASERECALIBRATOR } from '../modules/nf-core/gatk4/baserecalibrator/main'                                                                                                                                               
include { GATK4_APPLYBQSR } from '../modules/nf-core/gatk4/applybqsr/main'                                                                                                                                                             
include { GATK4_MUTECT2 } from '../modules/nf-core/gatk4/mutect2/main'                                                                                                                                                                 

workflow VARIANT_CALLING {
    take:
    bam
    bai
    fasta
    fai
    dict

    GATK4_SPLITNCIGARREADS(
        bam
            .filter{ meta, bam ->
                meta.bait != ""
            }.combine(bai, by:[0]),
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

}