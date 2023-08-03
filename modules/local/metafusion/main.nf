process METAFUSION {
    tag "$meta.id"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/metafusion:0.0.4' :
        'cmopipeline/metafusion:0.0.4' }"

    input:
    tuple val(meta), path(cff)
    path genebed
    path info
    path fasta
    path blocklist
    val numtools

    output:
    tuple val(meta), path("*final*cluster")             , emit: cluster
    tuple val(meta), path("*.filtered.cff")             , emit: filtered
    tuple val(meta), path("cis-sage.cluster")           , emit: cis
    tuple val(meta), path("problematic_chromosomes.cff"), emit: problem_chrom
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def sample = "${meta.sample}"
    """
    Metafusion_forte.sh  --cff $cff --outdir .  --gene_bed $genebed --gene_info $info  --genome_fasta $fasta --recurrent_bedpe $blocklist --num_tools=$numtools

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Metafusion docker: metafusion:0.0.4
        Metafusion_forte.sh: 0.0.1
    END_VERSIONS
    """
}
