include { STAR_ALIGN as STAR_FOR_ARRIBA     } from '../../modules/nf-core/star/align/main'
include { ARRIBA                            } from '../../modules/nf-core/arriba/main'
include { STAR_ALIGN as STAR_FOR_STARFUSION } from '../../modules/nf-core/star/align/main'
include { STARFUSION                        } from '../../modules/local/starfusion/detect/main'
include { FUSIONCATCHER_DETECT              } from '../../modules/local/fusioncatcher/detect/main'
include { FUSIONREPORT                      } from '../../modules/local/fusionreport/run/main'
include { ONCOKB_FUSIONANNOTATOR            } from '../../modules/local/oncokb/fusionannotator/main'
include { CSVTK_CONCAT as CSV_TO_TSV        } from '../../modules/nf-core/csvtk/concat/main'
include { TO_CFF as ARRIBA_TO_CFF} from '../../modules/local/convert_to_cff/main'
include { TO_CFF as FUSIONCATCHER_TO_CFF} from '../../modules/local/convert_to_cff/main'
include { TO_CFF as STARFUSION_TO_CFF} from '../../modules/local/convert_to_cff/main'
include { CSVTK_CONCAT as MERGE_CFF } from '../../modules/nf-core/csvtk/concat/main'
include {METAFUSION} from '../../modules/local/metafusion/main'


workflow FUSION {

    take:
    reads
    star_index
    gtf
    starfusion_ref
    fusioncatcher_ref
    fusion_report_db

    main:
    ch_versions = Channel.empty()
    fasta = params.fasta
    genebed = params.genebed
    info = params.info
    blocklist = params.blocklist


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
    MERGE_CFF(ARRIBA_TO_CFF.out.cff
            .map{ meta, file -> [meta, file]}
            .mix( 
                FUSIONCATCHER_TO_CFF.out.cff
                 .map{ meta, file -> [meta, file]}
            ).mix(
                STARFUSION_TO_CFF.out.cff
                 .map{ meta, file -> [meta, file]}
            ).groupTuple(by:[0]),
        'tsv',
        'tsv')

    METAFUSION(
        MERGE_CFF.out.csv,
        genebed,
        info,
        fasta,
        blocklist,
        "2"
    )

    ONCOKB_FUSIONANNOTATOR(METAFUSION.out.cluster)
    ch_versions = ch_versions.mix(ONCOKB_FUSIONANNOTATOR.out.versions.first())

    emit:
    ch_versions
}
