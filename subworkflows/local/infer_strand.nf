include { SEQTK_SAMPLE                } from '../../modules/nf-core/seqtk/sample/main'
include { PREPROCESS_READS            } from './preprocess_reads'
include { PICARD_COLLECTRNASEQMETRICS } from '../../modules/nf-core/picard/collectrnaseqmetrics/main'
include { EXTRACTSTRAND               } from '../../modules/local/extractstrand/main'
include { STAR_ALIGN                  } from '../../modules/nf-core/star/align/main'
include { GROUP_READS                 } from './group_reads'

workflow INFER_STRAND {

    take:
    reads
    star_index
    gtf
    refflat
    fasta

    main:

    ch_versions = Channel.empty()

    reads_branch = reads
        .branch{meta, reads ->
            auto: meta.strandedness == "auto"
            other: true
        }

    GROUP_READS(reads_branch.auto)

    SEQTK_SAMPLE(
        GROUP_READS.out.grouped_reads
            .map{ meta, reads ->
                [ meta, meta.single_end ? [reads[0]] : [reads[0], reads[1]], 50000 ]
            }
    )
    ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions.first())

    PREPROCESS_READS(SEQTK_SAMPLE.out.reads)
    ch_versions = ch_versions.mix(PREPROCESS_READS.out.ch_versions.first())

    STAR_ALIGN(
        PREPROCESS_READS.out.reads_untrimmed,
        star_index,
        gtf,
        false,
        [],
        []
    )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    PICARD_COLLECTRNASEQMETRICS(
        STAR_ALIGN.out.bam,
        refflat,
        fasta,
        []
    )
    ch_versions = ch_versions.mix(PICARD_COLLECTRNASEQMETRICS.out.versions.first())

    EXTRACTSTRAND(PICARD_COLLECTRNASEQMETRICS.out.metrics)

    amended_reads = EXTRACTSTRAND.out.strand
        .map{meta, strand_txt ->
            [ meta["sample"], strand_txt ]
        }.join(
            reads.map{ meta, reads ->
                [ meta["sample"], meta, reads ]
            }, by:[0]
        ).map{ sample, strand_txt, meta, reads ->
            def new_meta = meta.clone()
            new_meta["input_strandedness"] = new_meta["strandedness"]
            new_meta["strandedness"] = strand_txt.text.split("\\t")[2]
            [new_meta, reads]
        }.mix( reads_branch.other )

    emit:
    reads       = amended_reads
    ch_versions = ch_versions

}
