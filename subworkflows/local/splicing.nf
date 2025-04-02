include { PORTCULLIS_FULL } from '../../modules/nf-core/portcullis/full/main'

workflow SPLICING {
    take:
        bam
        junction_bed
        fasta
    main:
        ch_versions = Channel.empty()
        PORTCULLIS_FULL(
            bam,
            junction_bed,
            fasta
        )
        ch_versions = ch_versions.mix(PORTCULLIS_FULL.out.versions.first())
    emit:
        portcullis_bed = PORTCULLIS_FULL.out.pass_junctions_bed
        portcullis_tab = PORTCULLIS_FULL.out.pass_junctions_tab
        ch_versions = ch_versions

}
