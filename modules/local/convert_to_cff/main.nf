process TO_CFF {
    tag "$meta.id"
    label "process_single"

/// must be using singularity 3.7+
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' : 
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"
    
    input:
    tuple val(meta), val(caller), path(fusions)

    output:
    tuple val(meta), path("*.cff"), emit: cff

    when: 
    task.ext.when == null || task.ext.when

    script:
    def sample = "${meta.sample}"
    """
    make_cff_from_forte.R $caller $fusions $sample ${sample}_${caller}.cff

    """
}