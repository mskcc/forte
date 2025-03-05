include { GBCMS            } from '../../modules/msk/gbcms/main'
include { COMBINE_FILLOUTS } from '../../modules/local/combine_fillout_maf/main'
include { MAF_REFORMAT     } from '../../modules/local/reformat_fillout_maf/main'

workflow FILLOUT {

    take:
        bam
        bai
        maf
        fasta
        fai

    main:
    ch_versions = Channel.empty()

    maf_ch = bam
        .combine(maf)
        .filter{ meta, bam, meta2, maf ->
            meta.sample == meta2.sample
        }.map{ meta, bam, meta2, maf ->
            [ meta, maf ]
        }

    MAF_REFORMAT(maf_ch)
    ch_versions = ch_versions.mix(MAF_REFORMAT.out.versions.first())

    GBCMS(
        bam
            .combine(bai,by:[0])
            .combine(MAF_REFORMAT.out.maf,by:[0])
            .map{ meta, bam, bai, variants ->
                [ meta, bam, bai, variants, "${variants.getBaseName()}.gbcms.maf"]
            },
        fasta,
        fai

    )
    ch_versions = ch_versions.mix(GBCMS.out.versions.first())

    COMBINE_FILLOUTS(
        GBCMS.out.variant_file
            .combine(maf_ch,by:[0])
    )
    ch_versions = ch_versions.mix(COMBINE_FILLOUTS.out.versions.first())

    emit:
    combined_fillout_maf = COMBINE_FILLOUTS.out.maf
    ch_versions          = ch_versions

}
