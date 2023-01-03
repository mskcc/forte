include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTHSMETRICS     } from '../../modules/nf-core/picard/collecthsmetrics/main'
include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'

workflow QC {

    take:
    bam
    refflat
    rrna_intervals
    fastp_json

    main:
    fasta = params.fasta
    ch_versions = Channel.empty()

    PICARD_COLLECTRNASEQMETRICS(
        bam,
        refflat,
        fasta,
        rrna_intervals
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions.first())

    /*
    PICARD_COLLECTHSMETRICS(
        bam_ch.filter{ meta, bam ->
            meta.target != "none"
        },
        fasta,
        fai,
        bait,
        target
    )
    */

    multiqc_ch = PICARD_COLLECTRNASEQMETRICS.out.metrics
        .mix(fastp_json)
        .groupTuple(by:[0])

    MULTIQC(
        multiqc_ch.map{meta,multiqc_files -> multiqc_files},
        [],
        [],
        []
    )

    emit:
    ch_versions

}
