ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { PICARD_COLLECTHSMETRICS     } from '../../modules/nf-core/picard/collecthsmetrics/main'
include { MULTIQC                     } from '../../modules/nf-core/multiqc/main'
include { BAM_RSEQC                   } from '../nf-core/bam_rseqc/main'

workflow QC {

    take:
    bam
    bai
    refflat
    rrna_intervals
    rseqc_bed
    fastp_json
    htseq_counts
    star_log_final

    main:
    fasta = params.fasta
    ch_versions = Channel.empty()

    BAM_RSEQC(
        bam.join(bai, by:[0]),
        rseqc_bed,
        params.rseqc_modules
    )
    ch_versions = ch_versions.mix(BAM_RSEQC.out.versions.first())

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
        .mix(star_log_final)
        .mix(htseq_counts)
        .mix(BAM_RSEQC.out.bamstat_txt)
        .mix(BAM_RSEQC.out.innerdistance_freq)
        .mix(BAM_RSEQC.out.inferexperiment_txt)
        .mix(BAM_RSEQC.out.junctionannotation_log)
        .mix(BAM_RSEQC.out.junctionsaturation_rscript)
        .mix(BAM_RSEQC.out.readdistribution_txt)
        .mix(BAM_RSEQC.out.readduplication_pos_xls)
        .mix(BAM_RSEQC.out.tin_txt)
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
