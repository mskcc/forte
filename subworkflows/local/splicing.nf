include { PORTCULLIS_FULL } from '../../modules/nf-core/portcullis/full/main'
include { ONCO_JUNCS      } from '../../modules/local/oncogenic_junctions/main'

workflow SPLICING {
    take:
        bam
        junction_bed
        fasta
        reportable_junctions

    main:
        ch_versions = Channel.empty()
        PORTCULLIS_FULL(
            bam,
            junction_bed,
            fasta
        )

        ONCO_JUNCS(
            PORTCULLIS_FULL.out.pass_junctions_tab,
            reportable_junctions
        )

        ch_versions = ch_versions.mix(PORTCULLIS_FULL.out.versions.first())
        ch_versions = ch_versions.mix(ONCO_JUNCS.out.versions.first())

    emit:
        portcullis_bed = PORTCULLIS_FULL.out.pass_junctions_bed
        portcullis_tab = PORTCULLIS_FULL.out.pass_junctions_tab
        ch_versions = ch_versions

}
