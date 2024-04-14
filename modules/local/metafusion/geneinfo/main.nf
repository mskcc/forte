process METAFUSION_GENEINFO {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-28e4eeae4055c6901c8218fe9e55d342d59035ef:0b59d552b4c1c8b2845385236f1053a7e1fb80c6-0' :
        'quay.io/biocontainers/mulled-v2-28e4eeae4055c6901c8218fe9e55d342d59035ef:0b59d552b4c1c8b2845385236f1053a7e1fb80c6-0' }"

    input:
    tuple val(meta), path(gtf)
    path(starfusion_ref)
    path(fusioncatcher_ref)

    output:
    tuple val(meta), path("gene.info"), emit: metafusion_gene_info
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    make_gene_info_for_forte.R \\
        --primary_gtf $gtf \\
        --fc_custom_bed_gene_names $fusioncatcher_ref/custom_genes.bed \\
        --star_fusion_ref $starfusion_ref/ref_annot.gtf \\
        --fusioncatcher_ref $fusioncatcher_ref/organism.gtf \\
        --outputDir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        make_gene_info_for_forte.R: 0.0.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch gene.info

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        make_gene_info_for_forte.R: 0.0.1
    END_VERSIONS
    """
}
