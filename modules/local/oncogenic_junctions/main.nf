process ONCO_JUNCS{
    tag "$meta.id"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/tidyverse:4.4.2' :
        'ghcr.io/rocker-org/tidyverse:4.4.2' }"

    input:
    tuple val(meta), path(portcullis)
    path reportable_junctions

    output:
    tuple val(meta), path("*_oncogenic_isoforms.txt")        , emit: oncogenic_isoforms
    tuple val(meta), path("*_oncogenic_isoforms_dropped.txt"), emit: dropped_isoforms
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    oncogenic_isoforms.R \\
        --portcullis ${portcullis} \\
        --junctions ${reportable_junctions} \\
        --sample ${prefix} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        oncogenic_junctions.R: 0.0.1
    END_VERSIONS
    """
}