ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTHSMETRICS     } from '../../modules/nf-core/picard/collecthsmetrics/main'
include {
    MULTIQC ;
    MULTIQC as MULTIQC_COLLECT
} from '../../modules/nf-core/multiqc/main'

workflow QC {

    take:
    bam
    bai
    multiqc_files
    refflat
    rrna_intervals
    fai
    dict
    baits


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

    multiqc_files = multiqc_files
        .mix(PICARD_COLLECTRNASEQMETRICS.out.metrics)
        .mix(PICARD_COLLECTHSMETRICS.out.metrics)
        .map{ meta, file ->
            [meta.subMap(['sample']),file]
        }

    MULTIQC(
        multiqc_files.groupTuple(by:[0]),
        ch_multiqc_config.collect().ifEmpty([]),
        [],
        []
    )

    MULTIQC_COLLECT(
        multiqc_files.map{meta, multiqc_files -> multiqc_files}.collect().map{[[:],it]},
        ch_multiqc_config.collect().ifEmpty([]),
        [],
        []
    )

    emit:
    ch_versions

}
