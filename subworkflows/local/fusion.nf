include { STAR_ALIGN as STAR_FOR_ARRIBA     } from '../../modules/nf-core/star/align/main'
include { ARRIBA                            } from '../../modules/nf-core/arriba/main'
include { STAR_ALIGN as STAR_FOR_STARFUSION } from '../../modules/nf-core/star/align/main'
include { STARFUSION                        } from '../../modules/local/starfusion/detect/main'
include { FUSIONCATCHER_DETECT              } from '../../modules/local/fusioncatcher/detect/main'
include { ONCOKB_FUSIONANNOTATOR            } from '../../modules/local/oncokb/fusionannotator/main'
include { AGFUSION_BATCH                    } from '../../modules/local/agfusion/batch/main'
include { TO_CFF as ARRIBA_TO_CFF           } from '../../modules/local/convert_to_cff/main'
include { TO_CFF as FUSIONCATCHER_TO_CFF    } from '../../modules/local/convert_to_cff/main'
include { TO_CFF as STARFUSION_TO_CFF       } from '../../modules/local/convert_to_cff/main'
include { CAT_CAT as MERGE_CFF              } from '../../modules/nf-core/cat/cat/main'
include { METAFUSION                        } from '../../modules/local/metafusion/main'
include { ADD_FLAG                          } from '../../modules/local/add_flags/main'


workflow FUSION {

    take:
    reads
    star_index
    gtf
    starfusion_ref
    fusioncatcher_ref
    fusion_report_db
    agfusion_db
    pyensembl_cache

    main:
    ch_versions = Channel.empty()
    fasta = params.fasta
    gene_bed = params.metafusion_gene_bed
    gene_info = params.metafusion_gene_info
    blocklist = params.metafusion_blocklist

    STAR_FOR_ARRIBA(
        reads,
        star_index,
        gtf,
        false,
        [],
        []
    )
    ch_versions = ch_versions.mix(STAR_FOR_ARRIBA.out.versions.first())

    ARRIBA(
        STAR_FOR_ARRIBA.out.bam,
        fasta,
        gtf,
        [],
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(ARRIBA.out.versions.first())

    STAR_FOR_STARFUSION(
        reads,
        // use the star index in the starfusion reference to ensure compatibility
        starfusion_ref.map{ file( it + "/ref_genome.fa.star.idx")},
        gtf,
        false,
        [],
        []
    )
    ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions.first())

    //reads_junction = reads.join( STAR_FOR_STARFUSION.out.junction,by:[0] )

    STARFUSION(
        STAR_FOR_STARFUSION.out.junction.map{ meta, junction -> [ meta, [], junction ] },
        starfusion_ref
    )
    ch_versions = ch_versions.mix(STARFUSION.out.versions.first())

    FUSIONCATCHER_DETECT(
        reads,
        fusioncatcher_ref
    )
    ch_versions = ch_versions.mix(FUSIONCATCHER_DETECT.out.versions.first())

    fc_fusions = ["GRCh37","hg19","smallGRCh37"].contains(params.genome) ? FUSIONCATCHER_DETECT.out.fusions_alt : FUSIONCATCHER_DETECT.out.fusions


    ARRIBA_TO_CFF(ARRIBA.out.fusions
            .map{ meta, file ->[ meta, "arriba", file ] })
    FUSIONCATCHER_TO_CFF(fc_fusions
                    .map{ meta, file -> [ meta, "fusioncatcher", file ] } )
    STARFUSION_TO_CFF(STARFUSION.out.abridged
                    .map{ meta, file -> [ meta, "starfusion", file ] })
    // get expected number of callers for groupTuple
    numcallers = 1 + ( params.starfusion_url ? 1 : 0 ) + ( ["GRCh37","GRCh38"].contains(params.genome) ? 1 : 0 )

    MERGE_CFF(
        ARRIBA_TO_CFF.out.cff
            .map{ meta, file -> [meta, file]}
            .mix(
                FUSIONCATCHER_TO_CFF.out.cff
                    .map{ meta, file -> [meta, file]}
            ).mix(
                STARFUSION_TO_CFF.out.cff
                    .map{ meta, file -> [meta, file]}
            ).groupTuple(by:[0],size:numcallers),
    )

    METAFUSION(
        MERGE_CFF.out.file_out,
        gene_bed,
        gene_info,
        fasta,
        blocklist
    )

    ADD_FLAG(
        METAFUSION.out.cluster,
        METAFUSION.out.cis,
        METAFUSION.out.filtered,
        METAFUSION.out.problem_chrom
    )

    ONCOKB_FUSIONANNOTATOR(ADD_FLAG.out.cff)
    ch_versions = ch_versions.mix(ONCOKB_FUSIONANNOTATOR.out.versions.first())

    AGFUSION_BATCH(
        ADD_FLAG.out.cff,
        agfusion_db,
        pyensembl_cache
    )
    ch_versions = ch_versions.mix(AGFUSION_BATCH.out.versions.first())
    ch_versions = ch_versions.mix(ADD_FLAG.out.versions.first())
    ch_versions = ch_versions.mix(METAFUSION.out.versions.first())
    ch_versions = ch_versions.mix(ARRIBA_TO_CFF.out.versions.first())
    ch_versions = ch_versions.mix(FUSIONCATCHER_TO_CFF.out.versions.first())
    ch_versions = ch_versions.mix(STARFUSION_TO_CFF.out.versions.first())

    emit:
    ch_versions
}
