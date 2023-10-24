process COMBINE_FILLOUTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-tidyverse:2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/rocker-org/devcontainer/tidyverse:4' :
            'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta), path(fillout_maf), path(original_maf)

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo \$PATH
    rna_fillout_combine.R \\
        --fillout_maf ${fillout_maf}\\
        --maf ${original_maf} \\
        --output_maf ${original_maf.getBaseName()}.fillout.maf\\
        --rna_sample_id ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        rna_fillout_combine.R: 0.0.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.combined.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        rna_fillout_combine.R: 0.0.1
    END_VERSIONS
    """
}
