include { STAR_ALIGN as STAR_FOR_ARRIBA     } from '../../modules/nf-core/star/align/main'
include { ARRIBA                            } from '../../modules/nf-core/arriba/main'
include { STAR_ALIGN as STAR_FOR_STARFUSION } from '../../modules/nf-core/star/align/main'
include { STARFUSION                        } from '../../modules/local/starfusion/detect/main'
include { FUSIONCATCHER_DETECT              } from '../../modules/local/fusioncatcher/detect/main'
include { FUSIONREPORT                      } from '../../modules/local/fusionreport/run/main'

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

    FUSIONREPORT(
        ARRIBA.out.fusions
            .join(STARFUSION.out.abridged,by:[0])
            .join(fc_fusions,by:[0])
            .map{ meta, arriba, starfusion, fusioncatcher ->
                [
                    meta,
                    ["arriba","starfusion","fusioncatcher"], // names of callers recognized by fusionreport
                    [33,33,34], // weights
                    [arriba, starfusion, fusioncatcher] // files
                ]
            },
        fusion_report_db
    )

    emit:
    ch_versions
}
