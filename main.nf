#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mskcc/forte
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/mskcc/forte
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta          = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf            = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.starfusion_url = WorkflowMain.getGenomeAttribute(params, 'starfusion_url')
params.refflat        = WorkflowMain.getGenomeAttribute(params, 'refflat')
params.baits          = WorkflowMain.getGenomeAttribute(params, 'baits')
params.cdna           = WorkflowMain.getGenomeAttribute(params, 'cdna')
params.arriba_blacklist       = WorkflowMain.getGenomeAttribute(params, 'arriba_blacklist')
params.arriba_known_fusions   = WorkflowMain.getGenomeAttribute(params, 'arriba_known_fusions')
params.arriba_protein_domains = WorkflowMain.getGenomeAttribute(params, 'arriba_protein_domains')
params.blocklist       = WorkflowMain.getGenomeAttribute(params, 'block')
params.genebed  = WorkflowMain.getGenomeAttribute(params, 'genebed')
params.info = WorkflowMain.getGenomeAttribute(params, 'info')

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FORTE } from './workflows/forte'

//
// WORKFLOW: Run main mskcc/forte analysis pipeline
//
workflow MSKCC_FORTE {
    FORTE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    MSKCC_FORTE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
