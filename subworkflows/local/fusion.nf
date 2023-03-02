include { STAR_ALIGN as STAR_FOR_ARRIBA     } from '../../modules/nf-core/star/align/main'
include { ARRIBA                            } from '../../modules/nf-core/arriba/main'
include { STAR_ALIGN as STAR_FOR_STARFUSION } from '../../modules/nf-core/star/align/main'
include { STARFUSION                        } from '../../modules/local/starfusion/detect/main'
include { FUSIONCATCHER_DETECT              } from '../../modules/local/fusioncatcher/detect/main'
include { FUSIONREPORT                      } from '../../modules/local/fusionreport/run/main'
include { ONCOKB_FUSIONANNOTATOR            } from '../../modules/local/oncokb/fusionannotator/main'
include { CSVTK_CONCAT as CSV_TO_TSV        } from '../../modules/nf-core/csvtk/concat/main'

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

    STAR_FOR_ARRIBA(
        reads,
        star_index,
        gtf,
        false,
        false,
        false
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
        false,
        false
    )
    ch_versions = ch_versions.mix(STAR_FOR_STARFUSION.out.versions.first())

    reads_junction = reads.join( STAR_FOR_STARFUSION.out.junction )

    STARFUSION( reads_junction, starfusion_ref)
    ch_versions = ch_versions.mix(STARFUSION.out.versions.first())

    FUSIONCATCHER_DETECT(
        reads,
        fusioncatcher_ref
    )
    ch_versions = ch_versions.mix(FUSIONCATCHER_DETECT.out.versions.first())

    fc_fusions = ["GRCh37","hg19","smallGRCh37"].contains(params.genome) ? FUSIONCATCHER_DETECT.out.fusions_alt : FUSIONCATCHER_DETECT.out.fusions

    // get expected number of callers for groupTuple
    numcallers = 1 + ( params.starfusion_url ? 1 : 0 ) + ( ["GRCh37","GRCh38"].contains(params.genome) ? 1 : 0 )

    FUSIONREPORT(
        ARRIBA.out.fusions
            .map{ meta, file ->[ meta, "arriba", file ] }
            .mix(
                fc_fusions
                    .map{ meta, file -> [ meta, "fusioncatcher", file ] }
            ).mix(
                STARFUSION.out.abridged
                    .map{ meta, file -> [ meta, "starfusion", file ] }
            ).groupTuple(by:[0],size:numcallers)
            .map{ meta, caller, file ->
                def avg_weight = caller.collect({(100/caller.size()).toInteger()})
                avg_weight[-1] = avg_weight[-1] + (100-avg_weight.sum())
                [ meta, caller, avg_weight, file ]
            },
        fusion_report_db
    )
    ch_versions = ch_versions.mix(FUSIONREPORT.out.versions.first())

    CSV_TO_TSV(
        FUSIONREPORT.out.fusionreport_csv
            .map{ meta, csv ->
                [meta, [csv]]
            },
        "csv",
        "tsv"
    )
    ch_versions = ch_versions.mix(CSV_TO_TSV.out.versions.first())

    ONCOKB_FUSIONANNOTATOR(CSV_TO_TSV.out.csv)
    ch_versions = ch_versions.mix(ONCOKB_FUSIONANNOTATOR.out.versions.first())

    emit:
    ch_versions
}
