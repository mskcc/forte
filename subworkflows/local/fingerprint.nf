include { GBCMS            } from '../../modules/msk/gbcms/main'
include { FINGERPRINTPARSE } from '../../modules/local/fingerprintparse/main'

workflow FINGERPRINT {

    take:
    bam
    bai
    fp_vcf
    fasta
    fai

    main:
    ch_versions = Channel.empty()

    GBCMS(
        bam
            .combine(bai,by:[0])
            .combine(Channel.of(fp_vcf))
            .map{ meta, bam, bai, fp_vcf ->
                [ meta, bam, bai, fp_vcf, "${meta.id}.fp.vcf" ]
            },
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(GBCMS.out.versions)

    FINGERPRINTPARSE(GBCMS.out.variant_file)
    ch_versions = ch_versions.mix(FINGERPRINTPARSE.out.versions)

    emit:
    ch_versions
    fingerprint = FINGERPRINTPARSE.out.fp_txt
}
