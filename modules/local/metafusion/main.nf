process METAFUSION {
    tag "$meta.id"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/metafusion:0.0.6' :
        'cmopipeline/metafusion:0.0.6' }"

    input:
    tuple val(meta), path(cff)
    path genebed
    path info
    path fasta
    path blocklist

    output:
    tuple val(meta), path("*final*cluster")             , emit: cluster
    tuple val(meta), path("*.exons")                    , emit: cff
    tuple val(meta), path("cis-sage.cluster")           , emit: cis
    tuple val(meta), path("problematic_chromosomes.cff"), emit: problem_chrom
    tuple val(meta), path("filters.txt")                , emit: filters
    // tuple val(meta), path("*")                       , emit: all
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ""
    def sample = "${meta.sample}"
    """
    Metafusion_forte.sh \\
        --cff $cff \\
        --outdir . \\
        --gene_bed $genebed \\
        --gene_info $info \\
        --genome_fasta $fasta \\
        --recurrent_bedpe $blocklist \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Metafusion docker: \$METAFUSION_TAG
        Metafusion_forte.sh: 0.0.1
    END_VERSIONS
    """
}
