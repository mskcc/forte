process FUSIONREPORT {
    tag "$meta.id"
    label 'process_low'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda 'bioconda::star=2.7.9a'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/fusion-report:0.0.1' :
        'cmopipeline/fusion-report:0.0.1' }"

    input:
    tuple val(meta), val(callers), val(weights), path(fusions)
    path(fusionreport_ref)

    output:
    path "versions.yml"                                , emit: versions
    tuple val(meta), path("*fusionreport.tsv")         , emit: fusion_list
    tuple val(meta), path("*fusionreport_filtered.tsv"), emit: fusion_list_filtered
    tuple val(meta), path("*.html")                    , emit: report
    tuple val(meta), path("*.csv") , optional:true     , emit: fusionreport_csv
    tuple val(meta), path("*.json"), optional:true     , emit: fusionreport_json

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def caller_inputs = (1..callers.size()).collect{"--${callers[it-1]} ${fusions[it-1]} --${callers[it-1]}_weight ${weights[it-1]} "}.join(" ")
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fusion_report run \\
        $prefix \\
        . \\
        $fusionreport_ref \\
        $caller_inputs \\
        $args \\

    mv fusion_list.tsv ${prefix}.fusionreport.tsv
    mv fusion_list_filtered.tsv ${prefix}.fusionreport_filtered.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
