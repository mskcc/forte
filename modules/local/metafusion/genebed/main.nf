process METAFUSION_GENEBED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/rocker-org/devcontainer/tidyverse:4':
        'ghcr.io/rocker-org/devcontainer/tidyverse:4' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.metafusion.gene.bed"), emit: metafusion_gene_bed
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    final_generate_v75_gene_bed.R \\
        $gff \\
        ${prefix}.metafusion.gene.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        final_generate_v75_gene_bed.R: 0.0.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.metafusion.gene.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        final_generate_v75_gene_bed.R: 0.0.1
    END_VERSIONS
    """
}
