ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTHSMETRICS     } from '../../modules/nf-core/picard/collecthsmetrics/main'
include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'

workflow QC {

    take:
    bam
    bai
    refflat
    rrna_intervals
    fai
    dict
    baits
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

    PICARD_COLLECTHSMETRICS(
        bam
            .filter{ meta, bam ->
                meta.bait != ""
            }.combine(bai, by:[0])
            .combine(baits)
            .filter{ meta, bam, bai, bait, bait_file, target_file ->
                meta.bait == bait
            }.map{ meta, bam, bai, bait, bait_file, target_file ->
                [meta, bam, bai, bait_file, target_file]
            },
        [[:],fasta],
        fai,
        dict.map{ dict -> [[:],dict]}
    )

    multiqc_ch = PICARD_COLLECTRNASEQMETRICS.out.metrics
        .mix(fastp_json)
        .mix(PICARD_COLLECTHSMETRICS.out.metrics)
        .map{meta, multiqc_files -> multiqc_files }
        .collect()

    MULTIQC(
        multiqc_ch,
        ch_multiqc_config.collect().ifEmpty([]),
        [],
        []
    )

    emit:
    ch_versions

}
