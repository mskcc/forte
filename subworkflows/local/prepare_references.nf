include { STAR_GENOMEGENERATE            } from '../../modules/nf-core/star/genomegenerate/main'
include { UCSC_GTFTOGENEPRED             } from '../../modules/nf-core/ucsc/gtftogenepred/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST        } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { PREPARE_RRNA                   } from '../../modules/local/prepare_rrna/main'
include { GUNZIP                         } from '../../modules/nf-core/gunzip/main'
include { STARFUSION_DOWNLOAD            } from '../../modules/local/starfusion/download/main'

workflow PREPARE_REFERENCES {

    main:
    ch_versions = Channel.empty()

    if (params.gtf.endsWith(".gz")){
        GUNZIP([[:],params.gtf])
        gtf = GUNZIP.out.gunzip.map{ it[1] }
    } else {
        gtf = params.gtf
    }

    if (params.star_index){
        if (file(params.star_index).exists()){
            star_index = params.star_index
        } else {
            STAR_GENOMEGENERATE(params.fasta,gtf)
            star_index = STAR_GENOMEGENERATE.out.index
        }
    } else {
        STAR_GENOMEGENERATE(params.fasta,gtf)
        star_index = STAR_GENOMEGENERATE.out.index
    }

    UCSC_GTFTOGENEPRED(Channel.of([[id:params.genome],params.gtf]))
    ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)

    GATK4_CREATESEQUENCEDICTIONARY(params.fasta)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    PREPARE_RRNA(params.gtf)

    GATK4_BEDTOINTERVALLIST(
        PREPARE_RRNA.out.rRNA_bed.map{ bed -> [[id:params.genome],bed]},
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

    if (params.starfusion_reference) {
        if (file(params.starfusion_reference).exists()){
            starfusion_ref = params.starfusion_reference
        } else {
            STARFUSION_DOWNLOAD(params.starfusion_url)
            starfusion_ref = STARFUSION_DOWNLOAD.out.reference
        }
    } else {
        STARFUSION_DOWNLOAD(params.starfusion_url)
        starfusion_ref = STARFUSION_DOWNLOAD.out.reference
    }

    emit:
    star_index         = star_index
    // Convert queue channel to value channel so it never gets poison pilled
    refflat            = UCSC_GTFTOGENEPRED.out.refflat.map{it[1]}.first()
    reference_dict     = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    rrna_bed           = PREPARE_RRNA.out.rRNA_bed
    // Convert queue channel to value channel so it never gets poison pilled
    rrna_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list.map{it[1]}.first()
    gtf                = gtf
    starfusion_ref     = STARFUSION_DOWNLOAD.out.reference
    ch_versions        = ch_versions

}
