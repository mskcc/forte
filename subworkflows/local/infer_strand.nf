include { SAMTOOLS_VIEW               } from '../../modules/nf-core/samtools/view/main'
include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { EXTRACTSTRAND               } from '../../modules/local/extractstrand/main'


workflow INFER_STRAND {
    take:
    bam
    bai
    refflat
    fasta

    main:
    ch_versions = Channel.empty()

    bam = bam.branch{ meta, bam ->
        small: file(bam).size() < 10000000
        other: true
    }

    SAMTOOLS_VIEW(
        bam.other.join(bai,by:[0]),
        [[:],[]],
        []
    )

    PICARD_COLLECTRNASEQMETRICS(
        SAMTOOLS_VIEW.out.bam
            .mix(bam.other),
        refflat,
        fasta,
        []
    )

    EXTRACTSTRAND(PICARD_COLLECTRNASEQMETRICS.out.metrics)

    inferred_strand = EXTRACTSTRAND.out.strand
        .map{ meta, txt ->
            [meta, txt.text.split("\\t")[2]]
        }

    emit:
    ch_versions
    inferred_strand_txt = EXTRACTSTRAND.out.strand
    inferred_strand

}

workflow AMEND_STRAND {

    take:
    old_ch
    inferred_strandedness

    main:

    old_ch = old_ch
        .branch{
            auto: it[0]["strandedness"] == "auto"
            other: true
        }

    amended_ch = inferred_strandedness
        .join(old_ch.auto,by:[0])
        .map{ it ->
            def new_meta = it[0].clone()
            new_meta["input_strandedness"] = it[0]["strandedness"]
            new_meta["strandedness"] = it[1]
            [new_meta] + it[(2..-1)]
        }

    amended_ch = amended_ch.mix(old_ch.other)

    emit:
    amended_ch

}
