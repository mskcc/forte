process FUSION_FILTER {
    tag "$meta.id"
    label "process_single"

/// must be using singularity 3.7+
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/tidyverse:4.4.2' :
        'ghcr.io/rocker-org/tidyverse:4.4.2' }"

    input:
    tuple val(meta), path(cff), path(starfusion), path(fusioncatcher), path(arriba)
    path clinical_genes

    output:
    tuple val(meta), path("*_filtered_fusions.tsv")   , emit: filtered_fusions
    tuple val(meta), path("*_cis_sage_fusions.tsv")   , emit: cis_sage_fusions
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def sample = "${meta.sample}"
    """
    fusion_filtering.R \\
        --cff ${cff} \\
        --starfusion ${starfusion} \\
        --fusioncatcher ${fusioncatcher} \\
        --arriba ${arriba} \\
        --clinical_genes ${clinical_genes} \\
        --out_prefix ${sample} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        fusion_filtering.R: 0.0.1
    END_VERSIONS
    """
}
