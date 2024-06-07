process MAF_REFORMAT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-tidyverse:2.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'ghcr.io/rocker-org/devcontainer/tidyverse:4' :
            'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path("*.reformat.maf"), emit: maf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample}"
    """
    maf_reformat.R \\
        --maf ${maf}  \\
        --output_maf ${prefix}.reformat.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        maf_reformat.R: 0.0.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reformat.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        maf_reformat.R: 0.0.1
    END_VERSIONS
    """
}
