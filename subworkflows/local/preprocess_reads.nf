include { FASTP            } from '../../modules/nf-core/fastp/main'
include { UMITOOLS_EXTRACT } from '../../modules/nf-core/umitools/extract/main'

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

    reads_untrimmed = UMITOOLS_EXTRACT.out.reads
        .mix(
            reads.filter{ meta, reads -> ! meta.has_umi }
        )

    // once FASTP can run without producing reads then fix this logic.
    FASTP(
        reads_untrimmed,
        adapter_fasta,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    emit:
    reads_trimmed   = FASTP.out.reads
    reads_untrimmed = reads_untrimmed
    ungrouped_reads = trimmed_reads
    fastp_json      = FASTP.out.json
    ch_versions     = ch_versions
}
