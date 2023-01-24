process FUSIONREPORT {
    tag "$meta.id"
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/fusion-report:0.0.1' :
        'cmopipeline/fusion-report:0.0.1' }"

        //'docker.io/rannickscilifelab/fusion-report:2.1.5updated' :
        //'docker.io/rannickscilifelab/fusion-report:2.1.5updated' }"

    input:
    tuple val(meta), path(arriba_fusions), path(starfusion_fusions),  path(fusioncatcher_fusions)
    path(fusionreport_ref)

    output:
    path "versions.yml"                                , emit: versions
    tuple val(meta), path("*fusionreport.tsv")         , emit: fusion_list
    tuple val(meta), path("*fusionreport_filtered.tsv"), emit: fusion_list_filtered
    tuple val(meta), path("*.html")                    , emit: report
    tuple val(meta), path("*.csv"), optional:true      , emit: fusionreport_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def tools   = arriba_fusions ? "--arriba ${arriba_fusions} " : ''
    tools      += starfusion_fusions ? "--starfusion ${starfusion_fusions} " : ''
    tools      += fusioncatcher_fusions ? "--fusioncatcher ${fusioncatcher_fusions} " : ''
    def weights = "--starfusion_weight 33 --arriba_weight 33 --fusioncatcher_weight 34"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fusion_report run \\
        $meta.id \\
        . \\
        $fusionreport_ref \\
        $tools \\
        $weights \\
        --allow-multiple-gene-symbols \\
        --export csv

    mv fusion_list.tsv ${prefix}.fusionreport.tsv
    mv fusion_list_filtered.tsv ${prefix}.fusionreport_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
