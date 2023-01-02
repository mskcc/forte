#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_REFERENCES           } from '../../subworkflows/local/prepare_references'
include { TRIM_ALIGN as TRIM_ALIGN_UMI } from '../../subworkflows/local/trim_align'

params.run_alignment = true
params.run_qc = false
params.run_fusions = false

workflow test_trim_align {

    fq1   = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_1.fastq.gz"
    fq2   = "https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/fastq/test_rnaseq_2.fastq.gz"

    input_reads = [
        [id:'test', umi:'NNNXX', has_umi:true, single_end:false, strand:'reverse'],
        [ file(fq1), file(fq2) ]
    ]

    PREPARE_REFERENCES()

    TRIM_ALIGN_UMI(
        input_reads,
	PREPARE_REFERENCES.out.star_index,
	PREPARE_REFERENCES.out.gtf,
	true
    )
    
}
