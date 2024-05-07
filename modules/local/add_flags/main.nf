process ADD_FLAG {
    tag "$meta.id"
    label "process_single"

/// must be using singularity 3.7+
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' :
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta), path(cluster), path(cis), path(cff), path(problem_chrom), path(filters)
    path(clinicalgenes)

    output:
    tuple val(meta), path("*_metafusion_cluster.unfiltered.cff"), emit: unfiltered_cff
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample = "${meta.sample}"
    """
    add_flags_and_cluster_information.R \\
        $cff \\
        $cluster \\
        $cis \\
        $problem_chrom \\
        $filters \\
        $sample \\
        $clinicalgenes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        add_flags_and_cluster_information.R: 0.0.1
    END_VERSIONS
    """
}
