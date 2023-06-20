

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { TO_CFF as ARRIBA_TO_CFF} from './modules/local/convert_to_cff/main'
include { TO_CFF as FUSIONCATCHER_TO_CFF} from './modules/local/convert_to_cff/main'
include { TO_CFF as STARFUSION_TO_CFF} from './modules/local/convert_to_cff/main'
include { CSVTK_CONCAT as MERGE_CFF } from './modules/nf-core/csvtk/concat/main'
include {METAFUSION} from './modules/local/metafusion/main'
include { ONCOKB_FUSIONANNOTATOR            } from './modules/local/oncokb/fusionannotator/main'

nextflow.enable.dsl = 2

params.fasta          = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.blocklist       = WorkflowMain.getGenomeAttribute(params, 'block')
params.genebed  = WorkflowMain.getGenomeAttribute(params, 'genebed')
params.info = WorkflowMain.getGenomeAttribute(params, 'info')

WorkflowMain.initialise(workflow, params, log)
workflow {
    def meta1 = [:]
    meta1.id = "TESTING"
    meta1.sample  = "test"
    out1 = [meta1, "arriba", '/juno/work/ccs/noronhaa/tools/test_nf-core/anoronh4-forte_tests/full_test/results_review3/analysis/Sample_M17-9614_IGO_05500_GJ_12/arriba/Sample_M17-9614_IGO_05500_GJ_12.fusions.tsv']
    out2 = [meta1, "fusioncatcher",'/juno/work/ccs/noronhaa/tools/test_nf-core/anoronh4-forte_tests/full_test/results_review3/analysis/Sample_M17-9614_IGO_05500_GJ_12/fusioncatcher/Sample_M17-9614_IGO_05500_GJ_12.fusioncatcher.fusion-genes.hg19.txt']
    out3 = [meta1, "starfusion",'/juno/work/ccs/noronhaa/tools/test_nf-core/anoronh4-forte_tests/full_test/results_review3/analysis/Sample_M17-9614_IGO_05500_GJ_12/starfusion/Sample_M17-9614_IGO_05500_GJ_12.starfusion.fusion_predictions.tsv']
    ARRIBA_TO_CFF(out1)
    FUSIONCATCHER_TO_CFF(out2)
    STARFUSION_TO_CFF(out3)
    concat = ARRIBA_TO_CFF.out.cff
            .map{ meta, file -> [meta, file]}
            .mix( 
                FUSIONCATCHER_TO_CFF.out.cff
                 .map{ meta, file -> [meta, file]}
            ).mix(
                STARFUSION_TO_CFF.out.cff
                 .map{ meta, file -> [meta, file]}
            ).groupTuple(by:[0])

    MERGE_CFF(concat,
        'tsv',
        'tsv')

    METAFUSION(
        MERGE_CFF.out.csv,
        params.genebed,
        params.info,
        params.fasta,
        params.blocklist,
        "2"
    )

    ONCOKB_FUSIONANNOTATOR(METAFUSION.out.cluster)


    //TOCFF( ARRIBA.out.fusions
    ///      .map{ meta, file ->[ meta, "arriba", file ] })
}
