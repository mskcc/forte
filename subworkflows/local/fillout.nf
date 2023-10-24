
include { GETBASECOUNTSMULTISAMPLE } from '../../modules/local/getbasecountsmultisample/main'
include { COMBINE_FILLOUTS         } from '../../modules/local/combine_fillout_maf/main'

workflow FILLOUT {

    take:
        bam
        bai
        maf
        fasta
        fai

    main:
    ch_versions = Channel.empty()
    gbcms_ch = bam
        .combine(bai,by:[0])
        .combine(maf)
        .filter{ meta, bam, bai, meta2, maf ->
            meta.sample == meta2.sample
        }.map{ meta, bam, bai, meta2, maf ->
            [ meta, bam, bai, maf ]
        }

    GETBASECOUNTSMULTISAMPLE(gbcms_ch, fasta, fai)
    ch_versions = ch_versions.mix(GETBASECOUNTSMULTISAMPLE.out.versions.first())

    COMBINE_FILLOUTS(
        GETBASECOUNTSMULTISAMPLE.out.maf
            .combine(gbcms_ch, by:0)
            .view()
            .map{meta,fillout_maf,bam,bai,maf ->
                [meta,fillout_maf,maf]
            }
    )
    ch_versions = ch_versions.mix(COMBINE_FILLOUTS.out.versions.first())

    emit:
    combined_fillout_maf = COMBINE_FILLOUTS.out.maf
    ch_versions          = ch_versions

}
