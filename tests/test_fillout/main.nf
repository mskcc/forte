include { samplesheetToList } from 'plugin/nf-schema'
include { FILLOUT           } from '../../subworkflows/local/fillout'

workflow test_rna_fillout {

    // chr22 bam
    input_bam = [
        ['sample':'SAMPLE_PAIRED_END','id':'SAMPLE_PAIRED_END'],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam")
    ]
    input_bai = [
        ['sample':'SAMPLE_PAIRED_END','id':'SAMPLE_PAIRED_END'],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam.bai")
    ]

    // test maf has only chr22 variants
    fasta = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta")
    fai   = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta.fai")

    // run fillouts
    FILLOUT(
        Channel.of(input_bam),
        Channel.of(input_bai),
        Channel
            .fromList(samplesheetToList(params.maf_input, "assets/schema_maf_input.json")),
        fasta,
        fai
    )
}
