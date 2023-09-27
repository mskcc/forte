
include { MAF_INPUT_CHECK } from '../../subworkflows/local/input_check'
include { FILLOUT         } from '../../subworkflows/local/fillout'

workflow test_rna_fillout {

    // chr22 bam
    input_bam = [
        ['sample':'test','id':'test'],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam")
    ]
    input_bai = [
        ['sample':'test','id':'test'],
        file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.rna.paired_end.sorted.bam.bai")
    ]

    // test maf has only chr22 variants
    input_maf_samplesheet = "tests/test_fillout/data/input.tsv"
    fasta = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta")
    fai   = file("https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta.fai")

    // extract from input_maf_samplesheet
    MAF_INPUT_CHECK( input_maf_samplesheet )

    // run fillouts
    FILLOUT(
        Channel.of(input_bam),
        Channel.of(input_bai),
        MAF_INPUT_CHECK.out.mafs,
        fasta,
        fai
    )

}
