include { STAR_GENOMEGENERATE            } from '../../modules/nf-core/star/genomegenerate/main'
include { UCSC_GTFTOGENEPRED             } from '../../modules/nf-core/ucsc/gtftogenepred/main'
include { UCSC_GENEPREDTOBED             } from '../../modules/local/ucsc/genepredtobed/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST        } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { PREPARE_RRNA                   } from '../../modules/local/prepare_rrna/main'
include { GUNZIP                         } from '../../modules/nf-core/gunzip/main'
include { STARFUSION_DOWNLOAD            } from '../../modules/local/starfusion/download/main'
include { FUSIONCATCHER_DOWNLOAD         } from '../../modules/local/fusioncatcher/download/main'
include { FUSIONREPORT_DOWNLOAD          } from '../../modules/local/fusionreport/download/main'

workflow PREPARE_REFERENCES {

    main:
    ch_versions = Channel.empty()

    if (params.gtf.endsWith(".gz")){
        GUNZIP([[:],params.gtf])
        gtf = GUNZIP.out.gunzip.map{ it[1] }.first()
    } else {
        gtf = params.gtf
    }

    STAR_GENOMEGENERATE(params.fasta,gtf)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    star_index = STAR_GENOMEGENERATE.out.index

    UCSC_GTFTOGENEPRED(Channel.value(gtf).map{[[id:params.genome],it]})
    ch_versions = ch_versions.mix(UCSC_GTFTOGENEPRED.out.versions)

    UCSC_GENEPREDTOBED(UCSC_GTFTOGENEPRED.out.genepred)
    ch_versions = ch_versions.mix(UCSC_GENEPREDTOBED.out.versions)

    if ( params.refflat) {
        refflat = params.refflat
    } else {
        refflat = UCSC_GTFTOGENEPRED.out.refflat.map{it[1]}.first()
    }
    PREPARE_RRNA([],refflat)

    GATK4_CREATESEQUENCEDICTIONARY(params.fasta)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    GATK4_BEDTOINTERVALLIST(
        PREPARE_RRNA.out.rRNA_bed.map{ bed -> [[id:params.genome],bed]},
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

    if (params.starfusion_url) {
        STARFUSION_DOWNLOAD(params.starfusion_url)
        starfusion_ref = STARFUSION_DOWNLOAD.out.reference
    } else {
        starfusion_ref = Channel.empty()
    }

    if (["GRCh37","GRCh38"].contains(params.genome)){
        FUSIONCATCHER_DOWNLOAD()
        ch_versions = ch_versions.mix(FUSIONCATCHER_DOWNLOAD.out.versions)

        fusioncatcher_ref = FUSIONCATCHER_DOWNLOAD.out.reference
    } else {
        fusioncatcher_ref = Channel.empty()
    }

    cosmic_usr = params.cosmic_usr ?: ""
    cosmic_passwd = params.cosmic_passwd ?: ""
    FUSIONREPORT_DOWNLOAD(cosmic_usr,cosmic_passwd)


    emit:
    star_index         = star_index
    // Convert queue channel to value channel so it never gets poison pilled
    refflat            = refflat
    reference_dict     = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    rrna_bed           = PREPARE_RRNA.out.rRNA_bed
    // Convert queue channel to value channel so it never gets poison pilled
    rrna_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list.map{it[1]}.first()
    gtf                = gtf
    starfusion_ref     = starfusion_ref
    fusioncatcher_ref  = fusioncatcher_ref
    fusion_report_db   = FUSIONREPORT_DOWNLOAD.out.reference
    rseqc_bed          = UCSC_GENEPREDTOBED.out.bed.map{it[1]}.first()
    ch_versions        = ch_versions

}
