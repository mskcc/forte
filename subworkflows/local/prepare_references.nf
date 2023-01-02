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

    if (params.run_alignment) {
        STAR_GENOMEGENERATE(params.fasta,gtf)
        star_index = STAR_GENOMEGENERATE.out.index
    } else {
        star_index = Channel.ifEmpty([])
    }

    if (params.run_qc) {

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

	refflat = UCSC_GTFTOGENEPRED.out.refflat
	reference_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
	rrna_bed = PREPARE_RRNA.out.rRNA_bed
	rrna_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
    } else {
        refflat = Channel.empty()
	reference_dict = Channel.empty()
	rrna_bed = Channel.empty()
	rrna_interval_list = Channel.empty()

    }

    if (params.run_fusion){
        starfusion_ref = STARFUSION_DOWNLOAD(params.starfusion_reference)
    } else {
        starfusion_ref = Channel.empty()
    }

    emit:
    star_index         = star_index
    refflat            = refflat
    reference_dict     = reference_dict
    rrna_bed           = rrna_bed
    rrna_interval_list = rrna_interval_list
    gtf                = gtf
    starfusion_ref     = starfusion_ref
    ch_versions        = ch_versions

}
