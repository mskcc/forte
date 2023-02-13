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
                def rg_map = Utils.flowcellLaneFromFastq(meta.single_end ? reads : reads[0])
                meta_clone.read_group = "${meta.sample}@${rg_map["fcid"]}@${rg_map["lane"]}"
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
    if (params.skip_trimming){
        trimmed_reads = extracted_reads 
    } else {
        trimmed_reads = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }

    emit:
    reads           = trimmed_reads 
    fastp_json      = FASTP.out.json
    ch_versions     = ch_versions
}
