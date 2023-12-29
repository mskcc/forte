include { STAR_GENOMEGENERATE            } from '../../modules/nf-core/star/genomegenerate/main'
include { UCSC_GTFTOGENEPRED             } from '../../modules/nf-core/ucsc/gtftogenepred/main'
include { UCSC_GENEPREDTOBED             } from '../../modules/local/ucsc/genepredtobed/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { GATK4_BEDTOINTERVALLIST        } from '../../modules/nf-core/gatk4/bedtointervallist/main'
include { PREPARE_RRNA                   } from '../../modules/local/prepare_rrna/main'
include {
    GUNZIP as GUNZIP_GTF ;
    GUNZIP as GUNZIP_METAFUSIONGENEBED ;
    GUNZIP as GUNZIP_METAFUSIONBLOCKLIST
} from '../../modules/nf-core/gunzip/main'
include { STARFUSION_DOWNLOAD            } from '../../modules/local/starfusion/download/main'
include { FUSIONCATCHER_DOWNLOAD         } from '../../modules/local/fusioncatcher/download/main'
include { ARRIBA_DOWNLOAD                } from '../../modules/local/arriba/download/main'
include { KALLISTO_INDEX                 } from '../../modules/nf-core/kallisto/index/main'
include { AGFUSION_DOWNLOAD              } from '../../modules/local/agfusion/download/main'
include { AGAT_SPADDINTRONS              } from '../../modules/local/agat/spaddintrons/main'
include { METAFUSION_GENEBED             } from '../../modules/local/metafusion/genebed/main'
include { METAFUSION_GENEINFO            } from '../../modules/local/metafusion/geneinfo/main'

workflow PREPARE_REFERENCES {

    main:
    ch_versions = Channel.empty()

    if (params.gtf.endsWith(".gz")){
        GUNZIP_GTF([[:],params.gtf])
        gtf = GUNZIP_GTF.out.gunzip.map{ it[1] }.first()
    } else {
        gtf = params.gtf
    }

    if (params.metafusion_blocklist.endsWith(".gz")){
        GUNZIP_METAFUSIONBLOCKLIST([[:],params.metafusion_blocklist])
        metafusion_blocklist = GUNZIP_METAFUSIONBLOCKLIST.out.gunzip.map{ it[1] }.first()
    } else {
        metafusion_blocklist = params.metafusion_blocklist
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

    SAMTOOLS_FAIDX ([[:],params.fasta])

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

    ARRIBA_DOWNLOAD()

    AGAT_SPADDINTRONS(
        [[:],gtf],
        []
    )

    METAFUSION_GENEBED(
        AGAT_SPADDINTRONS.out.gff
    )

    METAFUSION_GENEINFO(
        [[:],gtf, starfusion_ref,fusioncatcher_ref]
    )

    AGFUSION_DOWNLOAD(
        params.ensembl_version,
        params.genome
    )
    ch_versions = ch_versions.mix(AGFUSION_DOWNLOAD.out.versions)

    KALLISTO_INDEX(params.cdna)
    ch_versions = ch_versions.mix(KALLISTO_INDEX.out.versions)

    emit:
    star_index         = star_index
    // Convert queue channel to value channel so it never gets poison pilled
    refflat            = refflat
    fasta_dict         = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    fasta_fai          = SAMTOOLS_FAIDX.out.fai
    rrna_bed           = PREPARE_RRNA.out.rRNA_bed
    // Convert queue channel to value channel so it never gets poison pilled
    rrna_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list.map{it[1]}.first()
    gtf                = gtf
    starfusion_ref     = starfusion_ref
    fusioncatcher_ref  = fusioncatcher_ref
    rseqc_bed          = UCSC_GENEPREDTOBED.out.bed.map{it[1]}.first()
    kallisto_index     = KALLISTO_INDEX.out.idx
    agfusion_db        = AGFUSION_DOWNLOAD.out.agfusion_db
    pyensembl_cache    = AGFUSION_DOWNLOAD.out.pyensembl_cache
    metafusion_blocklist = metafusion_blocklist
    metafusion_gene_bed = METAFUSION_GENEBED.out.metafusion_gene_bed
    metafusion_gene_info = METAFUSION_GENEINFO.out.metafusion_gene_info
    arriba_blacklist   = ARRIBA_DOWNLOAD.out.blacklist
    arriba_known_fusions = ARRIBA_DOWNLOAD.out.known_fusions
    arriba_protein_domains = ARRIBA_DOWNLOAD.out.protein_domains
    ch_versions        = ch_versions

}
