process CFF_ANNOTATE {
    tag "$meta.id"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' :
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta), path(cff), path(oncokb), path(agfusion)

    output:
    tuple val(meta), path("${prefix}.cff"), emit: cff
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${cff}"
    """
    add_annotations_cff.R \\
        --cff ${cff} \\
        --oncokb ${oncokb} \\
        --agfusion ${agfusion} \\
        --out-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        add_annotations_cff.R: 0.0.1
    END_VERSIONS
    """
}
