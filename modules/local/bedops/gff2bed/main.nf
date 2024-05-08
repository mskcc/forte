process BEDOPS_GFF2BED {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedops=2.4.41"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h9f5acd7_0' :
        'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    path('*.bed'), emit: bed
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gff2bed \\
       $args \\
       < $gff > ${prefix}.introns.bed
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedops: \$(bedops --version | awk '/version:/{print \$2}')
    END_VERSIONS
    """
}