process FUSIONREPORT_DOWNLOAD {
    tag 'fusionreport'
    label 'process_medium'

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.7.9a' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/fusion-report:0.0.1' :
        'cmopipeline/fusion-report:0.0.1' }"
        //'docker.io/rannickscilifelab/fusion-report:2.1.5updated' :
        //'docker.io/rannickscilifelab/fusion-report:2.1.5updated' }"

    output:
    path "db"                , emit: reference
    path "versions.yml"      , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    fusion_report download \\
        $args \\
        db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fusion_report: \$(fusion_report --version | sed 's/fusion-report //')
    END_VERSIONS
    """
}
