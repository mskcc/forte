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
    tuple val(meta), path("*final*cluster"), emit: cluster
    tuple val(meta), path("*.filtered.cff"), emit: filtered

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Metafusion_forte.sh  --cff $cff --outdir .  --gene_bed $genebed --gene_info  $info  --genome_fasta $fasta --recurrent_bedpe $blocklist --num_tools=$numtools

    """
}
