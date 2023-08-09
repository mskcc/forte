/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowForte.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.gtf, params.refflat, params.starfusion_url, params.metafusion_gene_bed, params.metafusion_blocklist, params.metafusion_gene_info ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/pipeline_multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { BAIT_INPUTS } from '../subworkflows/local/baits'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FORTE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: If baitsets are available, they will be added to the channel
    //
    BAIT_INPUTS ()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    PREPARE_REFERENCES()
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.ch_versions)


    PREPROCESS_READS(
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(PREPROCESS_READS.out.ch_versions)

    ALIGN_READS(
        PREPROCESS_READS.out.reads,
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
                PREPROCESS_READS.out.reads
                    .filter{ meta, reads -> ! ( meta.has_umi && params.dedup_umi_for_kallisto ) }
            ),
        PREPARE_REFERENCES.out.kallisto_index
    )
    ch_versions = ch_versions.mix(QUANTIFICATION.out.ch_versions)

    FUSION(
        PREPROCESS_READS.out.reads,
        PREPARE_REFERENCES.out.star_index,
        PREPARE_REFERENCES.out.gtf,
        PREPARE_REFERENCES.out.starfusion_ref,
        PREPARE_REFERENCES.out.fusioncatcher_ref,
        PREPARE_REFERENCES.out.fusion_report_db,
        PREPARE_REFERENCES.out.agfusion_db,
        PREPARE_REFERENCES.out.pyensembl_cache
    )
    ch_versions = ch_versions.mix(FUSION.out.ch_versions)

    QC_DEDUP(
        ALIGN_READS.out.bam_dedup,
        ALIGN_READS.out.bai_dedup,
        QUANTIFICATION.out.kallisto_log
            .filter{meta, log ->
                meta.has_umi && params.dedup_umi_for_kallisto
            },
        PREPARE_REFERENCES.out.refflat,
        PREPARE_REFERENCES.out.rrna_interval_list,
        PREPARE_REFERENCES.out.rseqc_bed,
        PREPARE_REFERENCES.out.fasta_fai,
        PREPARE_REFERENCES.out.fasta_dict,
        BAIT_INPUTS.out.baits
    )
    ch_versions = ch_versions.mix(QC_DEDUP.out.ch_versions)

    QC_DUP(
        ALIGN_READS.out.bam_withdup,
        ALIGN_READS.out.bai_withdup,
        PREPROCESS_READS.out.fastp_json
            .mix(QUANTIFICATION.out.htseq_counts)
            .mix(ALIGN_READS.out.star_log_final)
            .mix(
                QUANTIFICATION.out.kallisto_log
                    .filter{meta, log ->
                        ! (meta.has_umi && params.dedup_umi_for_kallisto)
                    }
            ),
        PREPARE_REFERENCES.out.refflat,
        PREPARE_REFERENCES.out.rrna_interval_list,
        PREPARE_REFERENCES.out.rseqc_bed,
        PREPARE_REFERENCES.out.fasta_fai,
        PREPARE_REFERENCES.out.fasta_dict,
        BAIT_INPUTS.out.baits
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowForte.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowForte.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    ch_multiqc_files = ch_multiqc_files.collect().map{ files -> [[:], files] }

    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
