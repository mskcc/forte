process ADD_FLAG {
    tag "$meta.id"
    label "process_single"

/// must be using singularity 3.7+
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/tidyverse:4.4.2' :
        'ghcr.io/rocker-org/tidyverse:4.4.2' }"

    input:
    tuple val(meta), path(cluster), path(cis), path(cff), path(problem_chrom), path(filters)
    path clinical_genes

    output:
    tuple val(meta), path("*_metafusion_cluster.unfiltered.cff")            , emit: unfiltered_cff
    tuple val(meta), path("*_metafusion_cluster.unfiltered.clinical.cff")   , emit: unfiltered_clinical_cff
    path "versions.yml"                                                     , emit: versions

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
        $sample

    cat *_metafusion_cluster.unfiltered.cff \\
        | head -1 > header.txt
    cat *_metafusion_cluster.unfiltered.cff \\
        | grep -iFwf $clinical_genes > tmp_clinicalgenes.txt
    cat header.txt tmp_clinicalgenes.txt > ${sample}_metafusion_cluster.unfiltered.clinical.cff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        add_flags_and_cluster_information.R: 0.0.1
    END_VERSIONS
    """
}
