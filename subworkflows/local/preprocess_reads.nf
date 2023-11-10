include { FASTP            } from '../../modules/nf-core/fastp/main'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/umitools/extract/main'
include {
    GROUP_READS as GROUP_READS_TRIMMED ;
    GROUP_READS as GROUP_READS_UNTRIMMED
} from './group_reads'

workflow PREPROCESS_READS {

    take:
    reads

    main:

    ch_versions = Channel.empty()

    adapter_fasta = file("${workflow.projectDir}/assets/adapters.fasta")
    if ( ! adapter_fasta.exists() ){ adapter_fasta = [] }

    reads = reads
        .map{ meta, reads ->
            def meta_clone = meta.clone()
            if (params.extract_fq_read_group) {
                def rg_map = Utils.flowcellLaneFromFastq(reads[0])
                meta_clone.read_group = "${meta.sample}@${rg_map["fcid"]}@${rg_map["lane"]}@${meta.fastq_pair_id}"
                meta_clone.id = meta_clone.read_group
            } else {
                meta_clone.read_group = meta_clone.id
            }
            [meta_clone, reads]
        }

    UMITOOLS_EXTRACT(
        reads.filter{ meta, reads -> meta.has_umi }
    )
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())

    // once FASTP can run without producing reads then fix this logic.
    FASTP(
        UMITOOLS_EXTRACT.out.reads
            .mix(
                reads.filter{ meta, reads -> ! meta.has_umi }
            ),
        adapter_fasta,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    GROUP_READS_UNTRIMMED(
        UMITOOLS_EXTRACT.out.reads
            .mix(
                reads.filter{ meta, reads -> ! meta.has_umi }
            )
    )
    reads_untrimmed = GROUP_READS_UNTRIMMED.out.grouped_reads

    GROUP_READS_TRIMMED(
        FASTP.out.reads
    )
    reads_trimmed = GROUP_READS_TRIMMED.out.grouped_reads

    emit:
    reads_trimmed   = reads_trimmed
    reads_untrimmed = reads_untrimmed
    fastp_json      = FASTP.out.json
    ch_versions     = ch_versions
}
