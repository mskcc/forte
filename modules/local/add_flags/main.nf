process ADD_FLAG{
    tag "$meta.id"
    label "process_single"

/// must be using singularity 3.7+
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' :
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta), path(cluster)
    tuple val(meta), path(cis)
    tuple val(meta), path(filtered)
    tuple val(meta), path(problem_chrom)

    output:
    tuple val(meta), path("*_metafusion_cluster.cff"), emit: cff
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample = "${meta.sample}"
    """
    add_flags_and_cluster_information.R $filtered $cluster $cis $problem_chrom $sample

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
    END_VERSIONS
    """
}