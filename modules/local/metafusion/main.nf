process METAFUSION {
    tag "$meta.id"
    label "process_low"

   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'cmopipeline/metafusion:0.0.2' : 
        'cmopipeline/metafusion:0.0.2' }"
    
    input:
    tuple val(meta), path(cff)
    path genebed
    path info
    path fasta
    path blocklist
    val numtools

    output:
    tuple val(meta), path("*cluster"), emit: cluster
    tuple val(meta), path("*.filtered.cff"), emit: filtered
    path "versions.yml"                          , emit: versions

    when: 
    task.ext.when == null || task.ext.when

    script:

    """
    Metafusion_forte.sh  --cff $cff --outdir .  --gene_bed $genebed --gene_info  $info  --genome_fasta $fasta --recurrent_bedpe $blocklist --num_tools=$numtools
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1)
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}
