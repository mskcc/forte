process METAFUSION_GENEBED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-28e4eeae4055c6901c8218fe9e55d342d59035ef:0b59d552b4c1c8b2845385236f1053a7e1fb80c6-0' :
        'quay.io/biocontainers/mulled-v2-28e4eeae4055c6901c8218fe9e55d342d59035ef:0b59d552b4c1c8b2845385236f1053a7e1fb80c6-0' }"

    input:
    tuple val(meta), path(gff)
    val ensembl_version

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
        ${ensembl_version}.metafusion.gene.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        final_generate_v75_gene_bed.R: 0.0.2
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
        final_generate_v75_gene_bed.R: 0.0.2
    END_VERSIONS
    """
}
