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

    extracted_reads = UMITOOLS_EXTRACT.out.reads
        .mix(
            reads.filter{ meta, reads -> ! meta.has_umi }
        )

    // once FASTP can run without producing reads then fix this logic.
    FASTP(
        extracted_reads,
        adapter_fasta,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    if (params.skip_trimming){
        trimmed_reads = extracted_reads
    } else {
        trimmed_reads = FASTP.out.reads
    }

    trimmed_grouped_reads = trimmed_reads
        .map{ meta, reads ->
            def read_group = meta.read_group
            def fastq_pair_id = meta.fastq_pair_id
            def meta_clone = meta.clone().findAll { !["read_group","fastq_pair_id"].contains(it.key) }
            meta_clone.id = meta.sample
            [groupKey(meta_clone,meta.fq_num), reads, read_group, fastq_pair_id]
        }.groupTuple(by:[0])
        .map{ meta, reads, read_group, fastq_pair_id ->
            meta = meta + [read_group:read_group.join(','), fastq_pair_id:fastq_pair_id.join(',')]
            [meta, reads.flatten()]
        }

    emit:
    reads           = trimmed_grouped_reads
    ungrouped_reads = trimmed_reads
    fastp_json      = FASTP.out.json
    ch_versions     = ch_versions
}
