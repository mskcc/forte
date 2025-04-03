/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAIT_INPUTS     } from '../subworkflows/local/baits'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PREPARE_REFERENCES                } from '../subworkflows/local/prepare_references'
include { PREPROCESS_READS                  } from '../subworkflows/local/preprocess_reads'
include { ALIGN_READS                       } from '../subworkflows/local/align_reads'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include {
    QC as QC_DUP ;
    QC as QC_DEDUP
} from '../subworkflows/local/qc'
include { EXTRACT_DEDUP_FQ                  } from '../subworkflows/local/extract_dedup_fq'
include { QUANTIFICATION                    } from '../subworkflows/local/quantification'
include { FUSION                            } from '../subworkflows/local/fusion'
include { FILLOUT                           } from '../subworkflows/local/fillout'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_forte_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FORTE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_maf_samplesheet // channel: samplesheet optionally read in from --maf_input
    main:

    ch_samplesheet = ch_samplesheet
        .groupTuple(by:[0])
        .map{ meta, reads ->
            def meta_clone = meta.clone()
            meta_clone.has_umi = meta.umi == [] ? false : true
            meta_clone.fq_num = reads.size()
            def fastq_pair_id = (1..reads.size()).toList().collect{ "${meta.id}_T${it}" }
            [meta_clone, reads, fastq_pair_id]
        }.transpose()
        .map{ meta, reads, fastq_pair_id ->
            def meta_clone = meta.clone()
            meta_clone.fastq_pair_id = fastq_pair_id
            [meta_clone,reads]
        }

    ch_versions = Channel.empty()

    ch_multiqc_files = Channel.empty()

    BAIT_INPUTS ()

    PREPARE_REFERENCES()
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.ch_versions)


    PREPROCESS_READS(
        ch_samplesheet
    )
    ch_versions = ch_versions.mix(PREPROCESS_READS.out.ch_versions)

    ALIGN_READS(
        params.skip_trimming ? PREPROCESS_READS.out.reads_untrimmed : PREPROCESS_READS.out.reads_trimmed,
        PREPARE_REFERENCES.out.star_index,
        PREPARE_REFERENCES.out.gtf
    )
    ch_versions = ch_versions.mix(ALIGN_READS.out.ch_versions)

    EXTRACT_DEDUP_FQ(
        ALIGN_READS.out.bam
            .filter{ meta, bam ->
                meta.has_umi && params.dedup_umi_for_kallisto
            }
    )
    ch_versions = ch_versions.mix(EXTRACT_DEDUP_FQ.out.ch_versions)

    QUANTIFICATION(
        ALIGN_READS.out.bam,
        ALIGN_READS.out.bai,
        PREPARE_REFERENCES.out.gtf,
        EXTRACT_DEDUP_FQ.out.dedup_reads
            .mix(
                params.skip_trimming ? PREPROCESS_READS.out.reads_untrimmed : PREPROCESS_READS.out.reads_trimmed
                    .filter{ meta, reads -> ! ( meta.has_umi && params.dedup_umi_for_kallisto ) }
            ),
        PREPARE_REFERENCES.out.kallisto_index
    )
    ch_versions = ch_versions.mix(QUANTIFICATION.out.ch_versions)

    FUSION(
        PREPROCESS_READS.out.reads_trimmed,
        PREPROCESS_READS.out.reads_untrimmed,
        PREPARE_REFERENCES.out.star_index,
        PREPARE_REFERENCES.out.fasta,
        PREPARE_REFERENCES.out.gtf,
        PREPARE_REFERENCES.out.starfusion_ref,
        PREPARE_REFERENCES.out.fusioncatcher_ref,
        PREPARE_REFERENCES.out.agfusion_db,
        PREPARE_REFERENCES.out.pyensembl_cache,
        PREPARE_REFERENCES.out.metafusion_gene_bed,
        PREPARE_REFERENCES.out.metafusion_gene_info,
        PREPARE_REFERENCES.out.metafusion_blocklist,
        workflow.profile.toString().split(",").contains("test") ? Channel.of([]).first() : PREPARE_REFERENCES.out.arriba_blacklist,
        workflow.profile.toString().split(",").contains("test") ? Channel.of([]).first() : PREPARE_REFERENCES.out.arriba_known_fusions,
        workflow.profile.toString().split(",").contains("test") ? Channel.of([]).first() : PREPARE_REFERENCES.out.arriba_protein_domains,
        params.clinical_genes,
        params.transcript_allowlist
    )
    ch_versions = ch_versions.mix(FUSION.out.ch_versions)

    FILLOUT(
        ALIGN_READS.out.bam,
        ALIGN_READS.out.bai,
        ch_maf_samplesheet,
        PREPARE_REFERENCES.out.fasta.map{ it[1] }.first(),
        PREPARE_REFERENCES.out.fasta_fai.map{ it[1] }.first()
    )
    ch_versions = ch_versions.mix(FILLOUT.out.ch_versions)

    QC_DEDUP(
        ALIGN_READS.out.bam_dedup,
        ALIGN_READS.out.bai_dedup,
        QUANTIFICATION.out.kallisto_log
            .mix(QUANTIFICATION.out.kallisto_count_feature)
            .filter{meta, log ->
                meta.has_umi && params.dedup_umi_for_kallisto
            }.mix(ALIGN_READS.out.umitools_dedup_log),
        PREPARE_REFERENCES.out.refflat,
        PREPARE_REFERENCES.out.rrna_interval_list,
        PREPARE_REFERENCES.out.rseqc_bed,
        PREPARE_REFERENCES.out.fasta,
        PREPARE_REFERENCES.out.fasta_fai,
        PREPARE_REFERENCES.out.fasta_dict,
        BAIT_INPUTS.out.baits
    )
    ch_versions = ch_versions.mix(QC_DEDUP.out.ch_versions)

    QC_DUP(
        ALIGN_READS.out.bam_withdup,
        ALIGN_READS.out.bai_withdup,
        PREPROCESS_READS.out.fastp_json
            .mix(ALIGN_READS.out.star_log_final)
            .mix(
                QUANTIFICATION.out.kallisto_log
                    .mix(QUANTIFICATION.out.kallisto_count_feature)
                    .filter{meta, log ->
                        ! (meta.has_umi && params.dedup_umi_for_kallisto)
                    }
            ),
        PREPARE_REFERENCES.out.refflat,
        PREPARE_REFERENCES.out.rrna_interval_list,
        PREPARE_REFERENCES.out.rseqc_bed,
        PREPARE_REFERENCES.out.fasta,
        PREPARE_REFERENCES.out.fasta_fai,
        PREPARE_REFERENCES.out.fasta_dict,
        BAIT_INPUTS.out.baits
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'forte_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
