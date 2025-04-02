include { GBCMS } from '../../modules/msk/gbcms/main'
include { FINGERPRINT_PARSE          } from '../../modules/local/fingerprint/parse/main'
include { FINGERPRINT_CORRELATION    } from '../../modules/local/fingerprint/correlation/main'

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

    FINGERPRINT_PARSE(GBCMS.out.variant_file)
    ch_versions = ch_versions.mix(FINGERPRINT_PARSE.out.versions)

    FINGERPRINT_CORRELATION(
        FINGERPRINT_PARSE.out.fp_txt
            .map{ meta, file -> file }
            .collect()
    )

    emit:
    ch_versions
    fingerprint = FINGERPRINT_PARSE.out.fp_txt
    fingerprint_correlation = FINGERPRINT_CORRELATION.out.correlation
}
